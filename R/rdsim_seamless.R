#' @title Risk Difference Simulation for Phase 2/3 Seamless Design
#' @description Simulate phase 2/3 seamless design testing for risk difference.
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
#' @param riskDiffH0s Numeric vector of length \eqn{M}. Risk differences
#'   under \eqn{H_0} for each active arm versus the common control. Defaults
#'   to 0 for superiority tests.
#' @param allocations Integer or integer vector of length \eqn{M + 1}.
#'   Number of subjects per arm within a randomization block. A single value
#'   implies equal allocation; defaults to 1.
#'   The first \eqn{M} elements refer to the active arms
#'   and the last element refers to the common control.
#' @param pis Numeric vector of length \eqn{M+1}. Each element corresponds
#'   to the response rate for a treatment arm. The first \eqn{M} elements
#'   refer to the active arms and the last element refers to the common control.
#' @param nullVariance 	Whether to use the variance under the null or the
#'   empirical variance under the alternative.
#' @param n Planned total sample size across all active arms and control.
#' @param plannedSubjects Numeric vector of length \eqn{K + 1} giving the planned
#'   cumulative sample size at each look.
#'   Each entry refers to the combined number of patients for the first
#'   active arm and the common control.
#' @param maxNumberOfIterations Number of Monte Carlo replications.
#'   Defaults to 1000.
#' @param seed Random seed for reproducibility.
#' @param nthreads Number of threads for parallel simulation.
#'   Use 0 to accept the default RcppParallel behavior.
#'
#' @return An S3 object of class \code{"rdsim_seamless"} with these components:
#'
#' * \code{overview}: A list summarizing trial-level results and settings:
#'     - \code{selectAsBest}: Probability of selecting each active arm as
#'       the best arm at the end of phase 2.
#'     - \code{selectToStage2}: Probability of selecting each active arm
#'       to enter stage 2.
#'     - \code{selectAnyToStage2}: Probability of selecting any active arm
#'       to enter stage 2.
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
#'     - \code{numberOfSubjects}: Cumulative enrollments by stage.
#'       replications that reached that stage.
#'     - \code{overallReject}: Overall probability of rejecting the null
#'       by trial end.
#'     - \code{overallFutility}: Overall probability of stopping for futility
#'       by trial end.
#'     - \code{expectedNumberOfEvents}: Expected cumulative events at trial end.
#'     - \code{expectedNumberOfSubjects}: Expected cumulative enrollments
#'       at trial end.
#'     - \code{criticalValues}: The input critical values for each stage.
#'     - \code{futilityBounds}: The input futility boundaries for each stage.
#'     - \code{riskDiffH0s}: The input risk differences under \eqn{H_0}.
#'     - \code{nullVariance}: Whether to use variance under \eqn{H_0}.
#'     - \code{numberOfIterations}: Number of simulation iterations performed.
#'     - \code{n}: Planned total sample size.
#'     - \code{allocations}: The input allocation ratios.
#'     - \code{responseRates}: The input response rates for each arm.
#'     - \code{plannedSubjects}: The input planned cumulative sample size at
#'       each look for the first active arm and the common control combined.
#'     - \code{M}: Number of active arms in Phase 2.
#'     - \code{K}: Number of sequential looks in Phase 3.
#'
#' * \code{sumdata1}: Data frame summarizing each iteration, stage, and
#'   treatment group:
#'     - \code{iterationNumber}, \code{stopStage}, \code{stageNumber},
#'       \code{treatmentGroup}, \code{accruals}, \code{events}, \code{phat}.
#'     - For each stage the final row summarizes the overall study
#'       (all arms combined).
#'
#' * \code{summdata2}: Data frame summarizing test statistics by iteration,
#'   stage, and active arm:
#'     - \code{iterationNumber}, \code{bestArm}, \code{stopStage},
#'       \code{stageNumber}, \code{activeArm},
#'       \code{totalAccruals}, \code{totalEvents},
#'       \code{riskDiff}, \code{vriskDiff}, \code{riskDiffZ},
#'       \code{reject}, \code{futility}.
#'     - For each active arm, total accruals and events refer to
#'       the combined counts for that arm and the common control at that stage.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' (sim1 <- rdsim_seamless(
#'   M = 3,
#'   K = 1,
#'   criticalValues = c(8, 2.349),
#'   futilityBounds = 1.036,
#'   pis = c(0.22, 0.25, 0.35, 0.20),
#'   n = 1120,
#'   plannedSubjects = c(280, 560),
#'   maxNumberOfIterations = 10000,
#'   seed = 314159,
#'   nthreads = 0))
#'
#' @export
rdsim_seamless <- function(
    M = 2,
    K = 1,
    criticalValues = NA,
    futilityBounds = NULL,
    riskDiffH0s = 0,
    allocations = 1,
    pis = NULL,
    nullVariance = TRUE,
    n = NA,
    plannedSubjects = NA,
    maxNumberOfIterations = 1000,
    seed = 0,
    nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  rdsim_seamless_Rcpp(
    M, K, criticalValues, futilityBounds, riskDiffH0s,
    allocations, pis, nullVariance, n, plannedSubjects,
    maxNumberOfIterations, seed)
}
