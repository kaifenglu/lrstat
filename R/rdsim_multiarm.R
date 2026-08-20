#' @title Risk Difference Simulation for Multi-Arm Multi-Stage Design
#' @description Simulate a multi-arm multi-stage design for binary responses
#'   using risk-difference Wald statistics with closed testing.
#'
#' @param M Number of active treatment arms.
#' @param kMax Number of sequential looks.
#' @param criticalValues Numeric matrix of dimension
#'   \eqn{kMax \times M} giving the by-look critical values for the closed
#'   testing procedure. The first column is used for the level-M test and the
#'   last column for the level-1 test.
#' @param futilityBounds Numeric vector of length \eqn{kMax - 1} giving the
#'   futility boundaries on the Wald-statistic scale for the first
#'   \eqn{kMax - 1} looks. At an interim look, the study stops for futility if
#'   all active treatment arms fall below the futility boundary. If omitted,
#'   no interim futility stopping is applied.
#' @param riskDiffH0s Scalar or numeric vector of length \eqn{M}. Risk differences
#'   under \eqn{H_0} for each active arm versus the common control. Defaults to 0.
#' @param allocations Integer or integer vector of length \eqn{M + 1}. Number
#'   of subjects per arm within a randomization block. A single value implies
#'   equal allocation; defaults to 1. The first \eqn{M} elements refer to the
#'   active arms and the last element refers to the common control.
#' @param pis Numeric vector of length \eqn{M + 1}. Each element corresponds to
#'   the response rate for a treatment arm. The first \eqn{M} elements refer to
#'   the active arms and the last element refers to the common control.
#' @param nullVariance Whether to use the variance under the null or the
#'   empirical variance under the alternative.
#' @param n Planned total sample size across all active arms and control.
#' @param plannedSubjects Numeric vector of length \eqn{kMax} giving the planned
#'   cumulative sample size at each look for the first active arm and the
#'   common control combined.
#' @param maxNumberOfIterations Number of Monte Carlo replications.
#'   Defaults to 1000.
#' @param seed Random seed for reproducibility.
#' @param nthreads Number of threads for parallel simulation.
#'   Use 0 to accept the default RcppParallel behavior.
#'
#' @return An S3 object of class \code{"rdsim_multiarm"} with these components:
#'
#' * \code{overview}: A list summarizing trial-level results and settings:
#'     - \code{overallReject}: Overall probability of rejecting the null
#'       by trial end.
#'     - \code{overallFutility}: Overall probability of stopping for futility
#'       by trial end.
#'     - \code{rejectPerStage}: Probability of rejecting the null for each
#'       active arm at each stage.
#'     - \code{futilityPerStage}: Probability of futility stopping for each
#'       active arm at each stage.
#'     - \code{cumulativeRejection}: Cumulative probability of rejection by
#'       stage.
#'     - \code{cumulativeFutility}: Cumulative futility stopping probability
#'       by stage.
#'     - \code{numberOfEvents}: Cumulative event counts by stage and arm.
#'     - \code{numberOfSubjects}: Cumulative enrollments by stage and arm.
#'     - \code{expectedNumberOfEvents}: Expected cumulative events at trial end.
#'     - \code{expectedNumberOfSubjects}: Expected cumulative enrollments at
#'       trial end.
#'     - \code{criticalValues}: The input matrix of by-level critical values.
#'     - \code{futilityBounds}: The input futility boundaries for each stage.
#'     - \code{riskDiffH0s}: The input risk differences under \eqn{H_0}.
#'     - \code{nullVariance}: Whether to use variance under \eqn{H_0}.
#'     - \code{numberOfIterations}: Number of simulation iterations performed.
#'     - \code{n}: Planned total sample size.
#'     - \code{allocations}: The input allocation ratios.
#'     - \code{responseRates}: The input response rates for each arm.
#'     - \code{plannedSubjects}: The input planned cumulative sample size at
#'       each look for the first active arm and the common control combined.
#'     - \code{M}: Number of active arms.
#'     - \code{kMax}: Number of sequential looks.
#'
#' * \code{sumdata1}: Data frame summarizing each iteration, stage, and
#'   treatment group:
#'     - \code{iterationNumber}, \code{stopStage}, \code{stageNumber},
#'       \code{treatmentGroup}, \code{accruals}, \code{events}, \code{phat}.
#'     - For each stage the final row summarizes the overall study
#'       (all arms combined).
#'
#' * \code{sumdata2}: Data frame summarizing test statistics by iteration,
#'   stage, and active arm:
#'     - \code{iterationNumber}, \code{stopStage}, \code{stageNumber},
#'       \code{activeArm}, \code{totalAccruals}, \code{totalEvents},
#'       \code{riskDiff}, \code{vriskDiff}, \code{riskDiffZ},
#'       \code{reject}, \code{futility}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' (sim1 <- rdsim_multiarm(
#'   M = 2,
#'   kMax = 3,
#'   criticalValues = matrix(c(3.880, 2.747, 2.275,
#'                             3.710, 2.511, 1.993), 3, 2),
#'   futilityBounds =  c(0.043, 1.194),
#'   pis = c(0.25, 0.30, 0.20),
#'   n = 486,
#'   plannedSubjects = c(146, 292, 324),
#'   maxNumberOfIterations = 10000,
#'   seed = 314159,
#'   nthreads = 0))
#'
#' @export
rdsim_multiarm <- function(
    M = 2,
    kMax = 1,
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

  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  rdsim_multiarm_Rcpp(
    M, kMax, criticalValues, futilityBounds, riskDiffH0s,
    allocations, pis, nullVariance, n, plannedSubjects,
    maxNumberOfIterations, seed)
}
