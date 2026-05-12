#' @title Log-Rank Test Simulation for Phase 2/3 Seamless Design
#' @description Simulate phase 2/3 seamless design using a weighted
#' log-rank test. Analyses can be triggered either by the cumulative number
#' of events (combined for an active arm and the common control) or by
#' pre-specified calendar times.
#'
#' @param M Number of active treatment arms in Phase 2.
#' @param K Number of sequential looks in Phase 3.
#' @param criticalValues Numeric vector of length \eqn{K + 1} giving the
#'   critical value for the Wald statistic at each look (Look 1 through
#'   Look \eqn{K + 1}). Decision rule:
#'   - At Look 1, compute the Wald statistic for each active arm versus the
#'     common control. If the maximum of these statistics exceeds the Look 1
#'     critical value, stop for efficacy.
#'   - If the Look 1 stopping rule is not met, select the active arm with the
#'     largest Wald statistic as the "best" arm and continue with that arm
#'     only versus control at subsequent looks.
#'   - For each look \eqn{j = 2,\ldots,K+1}, compare the selected arm to
#'     control; if its Wald statistic exceeds the Look \eqn{j} critical
#'     value, stop for efficacy; otherwise continue.
#'   - If no critical value is exceeded by Look \eqn{K + 1}, the procedure
#'     ends without rejection.
#' @param futilityBounds Numeric vector of length \eqn{K} giving the futility
#'   boundaries for Phase 2 and the first \eqn{K-1} looks in Phase 3.
#'   The study stops for futility:
#'   - in Phase 2 if all active treatment arms cross the phase-2 futility
#'     boundary;
#'   - in Phase 3 if the selected arm crosses the futility boundary at an
#'     interim look;
#'   If omitted, no interim futility stopping is applied.
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
#'     - \code{selectAsBest}: Probability of selecting each active arm as
#'       the best arm at the end of phase 2.
#'     - \code{rejectPerStage}: Probability of rejecting the null for each
#'       active arm at each stage.
#'     - \code{futilityPerStage}: Probability of futility stopping for each
#'       active arm at each stage.
#'     - \code{cumulativeRejection}: Cumulative probability of rejection by stage.
#'     - \code{cumulativeFutility}: Cumulative futility stopping probabilities
#'       by stage.
#'     - \code{numberOfEvents}: Cumulative event counts by stage, including
#'       events from all arms in stage 1 and events from the selected arm
#'       and control in later stages.
#'     - \code{numberOfDropouts}: Cumulative dropouts by stage.
#'     - \code{numberOfSubjects}: Cumulative enrollments by stage.
#'     - \code{analysisTime}: Average calendar time for each stage among
#'       replications that reached that stage.
#'     - \code{overallReject}: Overall probability of rejecting the null
#'       by trial end.
#'     - \code{overallFutility}: Overall probability of stopping for futility
#'       by trial end.
#'     - \code{expectedNumberOfEvents}: Expected cumulative events at trial end.
#'     - \code{expectedNumberOfDropouts}: Expected cumulative dropouts at trial end.
#'     - \code{expectedNumberOfSubjects}: Expected cumulative enrollments
#'       at trial end.
#'     - \code{expectedStudyDuration}: Expected study duration.
#'     - \code{criticalValues}: The input critical values for each stage.
#'     - \code{futilityBounds}: The input futility boundaries for each stage.
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
#'     - \code{iterationNumber}, \code{bestArm}, \code{stopStage},
#'       \code{stageNumber}, \code{analysisTime}, \code{activeArm},
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
#' (sim1 <- lrsim_seamless(
#'   M = 2,
#'   K = 2,
#'   criticalValues = c(3.882, 2.733, 2.222),
#'   futilityBounds = c(0.259, 1.201),
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   lambdas = list(log(2)/12*0.5, log(2)/12*0.7, log(2)/12),
#'   n = 700,
#'   plannedEvents = c(42, 84, 126),
#'   maxNumberOfIterations = 10000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 0))
#'
#' @export
lrsim_seamless <- function(
    M = 2,
    K = 1,
    criticalValues = NA,
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

  lrsim_seamless_Rcpp(
    M, K, criticalValues, futilityBounds, hazardRatioH0s,
    allocations, accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction, lambdas, gammas,
    n, followupTime, fixedFollowup, rho1, rho2,
    plannedEvents, plannedTime, maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed)
}
