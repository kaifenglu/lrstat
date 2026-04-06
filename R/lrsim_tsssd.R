#' @title Log-Rank Test Simulation for Multiple Active Arms vs. Common Control
#' @description Performs simulation for multiple-arm group sequential trials
#' based on weighted log-rank test. The looks are driven by the total
#' number of events in Arm 1 and common control combined. Alternatively,
#' the analyses can be planned to occur at specified calendar times.
#'
#' @param M The number of active treatment arms in phase 2. Must be at least 2.
#' @param K The number of looks in phase 3. Must be at least 1.
#' @param criticalValues A numeric vector of length \eqn{K + 1} that
#'   gives the critical values for the Wald statistic at each look
#'   (Look 1 through Look \eqn{K + 1}). Decision rule:
#'
#'   - At Look 1, compute the Wald statistic for each active arm versus the
#'     common control. If the maximum of these statistics exceeds the
#'     Look 1 critical value, the trial stops for efficacy.
#'
#'   - If the Look 1 stopping rule is not met, select the active arm with
#'     the largest Wald statistic as the best arm and proceed to Look 2
#'     (comparing only that selected arm to the common control).
#'
#'   - For each subsequent look \eqn{j = 2,\ldots, K + 1}, compare the
#'     selected active arm to the common control; if its Wald statistic
#'     exceeds the critical value for Look \eqn{j} the trial stops for
#'     efficacy, otherwise continue to the next look.
#'
#'   - If no critical value is exceeded by Look \code{K + 1}, the procedure
#'     ends without rejection.
#' @param hazardRatioH0s Hazard ratios under the null hypothesis for
#'   each active arm versus the common control. A vector of length M. Defaults
#'   to 1 for superiority test.
#' @param allocations Number of subjects in each arm in a randomization block.
#'   Defaults to 1 for equal randomization.
#' @inheritParams param_accrualTime
#' @inheritParams param_accrualIntensity
#' @inheritParams param_piecewiseSurvivalTime
#' @inheritParams param_stratumFraction
#' @param lambdas A list of vectors of hazard rates for event in each
#'   analysis time interval by stratum for each arm.
#' @param gammas A list of vectors of hazard rates for dropout in each
#'   analysis time interval by stratum for each arm.
#' @param n Sample size.
#' @inheritParams param_followupTime
#' @inheritParams param_fixedFollowup
#' @inheritParams param_rho1
#' @inheritParams param_rho2
#' @param plannedEvents The planned cumulative total number of events at
#'   Look 1 to Look \code{K + 1} for Arms A and C combined.
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
#' @return A list with 3 components:
#'
#' * \code{overview}: A list containing the following information:
#'
#'     - \code{selectAsBest}: The probability of selecting each active arm
#'       as the best arm at phase 2 analysis.
#'
#'     - \code{rejectPerStage}: The probability of rejecting the null
#'       hypothesis for each active arm at each stage.
#'
#'     - \code{cumulativeRejection}: The cumulative probability of rejecting
#'       the null hypothesis for each active arm by stage.
#'
#'     - \code{numberOfEvents}: The cumulative number of events up to each
#'       stage, including events from all arms in stage 1 and new
#'       events from the selected active arm and the common control
#'       in later stages.
#'
#'     - \code{numberOfDropouts}: The cumulative number of dropouts up to each
#'       stage, including dropouts from all arms in stage 1 and new
#'       dropouts from the selected active arm and the common control
#'       in later stages.
#'
#'     - \code{numberOfSubjects}: The cumulative number of subjects enrolled
#'       up to each stage, including subjects from all arms in stage 1 and new
#'       subjects from the selected active arm and the common control in
#'       later stages.
#'
#'     - \code{analysisTime}: The average calendar time for each stage among
#'       replications that have the stage (do not reject before the stage).
#'
#'     - \code{overallReject}: The overall probability of rejecting the
#'       null hypothesis by the end of the trial.
#'
#'     - \code{expectedNumberOfEvents}: The expected cumulative number of
#'       events at the end of the trial.
#'
#'     - \code{expectedNumberOfDropouts}: The expected cumulative number of
#'       dropouts at the end of the trial.
#'
#'     - \code{expectedNumberOfSubjects}: The expected cumulative number of
#'       subjects enrolled at the end of the trial.
#'
#'     - \code{expectedStudyDuration}: The expected study duration.
#'
#'     - \code{hazardRatioH0s}: The hazard ratios under the null hypothesis
#'       for each active arm versus the common control.
#'
#'     - \code{useEvents}: Whether the analyses are planned based on the
#'       number of events.
#'
#'     - \code{numberOfIterations}: The number of simulation iterations.
#'
#'     - \code{n}: The planned sample size for the trial.
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
#'     - \code{M}: The number of active arms in phase 2.
#'
#'     - \code{K}: The number of looks in phase 3.
#'
#' * \code{sumdata1}: A data frame of summary data by iteration, stage, and arm:
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
#'     - \code{treatmentGroup}: The treatment group. For each stage, the final
#'       row corresponds to the overall study (all arms combined).
#'
#'     - \code{accruals}: The number of subjects enrolled at the stage for
#'       the treatment group.
#'
#'     - \code{events}: The number of events at the stage for the
#'       treatment group.
#'
#'     - \code{dropouts}: The number of dropouts at the stage for the
#'       treatment group.
#'
#' * \code{summdata2}: A data frame of summary data by iteration, stage, and
#'   active arm for the log-rank test statistics comparing the active arm to
#'   the common control:
#'
#'     - \code{iterationNumber}: The iteration number.
#'
#'     - \code{stageNumber}: The stage number, covering all stages even if
#'       the trial stops at an interim look.
#'
#'     - \code{analysisTime}: The time for the stage since trial start.
#'
#'     - \code{activeArm}: The active arm under consideration.
#'
#'     - \code{totalAccruals}: The total number of subjects enrolled at
#'       the stage for the active arm and the common control combined.
#'
#'     - \code{totalEvents}: The total number of events at the stage for
#'       the active arm and the common control combined.
#'
#'     - \code{totalDropouts}: The total number of dropouts at the stage for
#'       the active arm and the common control combined.
#'
#'     - \code{uscore}: The log-rank test score statistic comparing the
#'       active arm to the control.
#'
#'     - \code{vscore}: The log-rank test variance statistic comparing the
#'       active arm to the control.
#'
#'     - \code{logRankStatistic}: The log-rank test Z-statistic comparing the
#'       active arm to the control.
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
#' sim1 = lrsim_tsssd(
#'   M = 2,
#'   K = 2,
#'   criticalValues = c(3.852050, 2.723811, 2.223982),
#'   accrualTime = c(0, 8),
#'   accrualIntensity = c(10, 28),
#'   piecewiseSurvivalTime = 0,
#'   lambdas = list(log(2)/12*0.5, log(2)/12*0.7, log(2)/12),
#'   n = 700,
#'   plannedEvents = c(40, 80, 120),
#'   maxNumberOfIterations = 10000,
#'   maxNumberOfRawDatasetsPerStage = 1,
#'   seed = 314159,
#'   nthreads = 0)
#'
#' head(sim1$sumdata1)
#' head(sim1$sumdata2)
#' head(sim1$rawdata)
#'
#' @export
lrsim_tsssd <- function(
    M = 2,
    K = 1,
    criticalValues = NA,
    hazardRatioH0s = 1,
    allocations = 1,
    accrualTime = 0,
    accrualIntensity = NA,
    piecewiseSurvivalTime = 0,
    stratumFraction = 1,
    lambdas = list(NA, NA, NA),
    gammas = list(0, 0, 0),
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

  lrsim_tsssd_Rcpp(
    M, K, criticalValues, hazardRatioH0s,
    allocations, accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction, lambdas, gammas,
    n, followupTime, fixedFollowup, rho1, rho2,
    plannedEvents, plannedTime, maxNumberOfIterations,
    maxNumberOfRawDatasetsPerStage, seed)
}
