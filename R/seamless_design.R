#' @title Exit Probabilities for Phase 2/3 Seamless Design
#' @description Computes the upper and lower exit probabilities for a phase
#' 2/3 seamless design. In Phase 2, multiple active arms are compared
#' against a common control arm. If the test statistic for the arm
#' ranked \code{rankp0} at the end of Phase 2 crosses the
#' efficacy boundary, the trial stops early for efficacy; if it falls below
#' the futility boundary, the trial stops early for futility. Otherwise,
#' the arm is selected to proceed to Phase 3, where it is tested against
#' the control over multiple looks with upper and optional
#' lower stopping boundaries.
#'
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param theta A vector of length \eqn{M} representing the true treatment
#'   effects for each active arm versus the common control.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#' @param K Number of sequential looks in Phase 3.
#' @param b A vector of efficacy boundaries (length \eqn{K + 1}). The first
#'   element is the efficacy boundary for the Phase-2 test statistic;
#'   the remaining \eqn{K} elements are efficacy boundaries for the selected
#'   arm in Phase 3.
#' @param a An optional vector of futility boundaries (length \eqn{K + 1}).
#'   The first element is the futility boundary for the Phase-2 test statistic;
#'   the remaining \eqn{K} elements are futility boundaries for the selected
#'   arm in Phase 3. If omitted, no futility stopping is applied.
#' @param I A vector of information levels (length \eqn{K + 1}) for any active
#'   arm versus the common control. The first element is for Phase 2;
#'   the remaining \eqn{K} elements are for the looks in Phase 3.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward when the trial continues to Phase 3.
#'   \code{rankp0 = 1} selects the largest Phase-2 Z-statistic,
#'   \code{rankp0 = 2} selects the second largest, and so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' The function assumes a multivariate normal distribution for the Wald
#' statistics. Among designs that continue beyond the Phase-2 analysis,
#' the carried-forward arm is the one with rank \code{rankp0} based on the
#' p-value of the Z-statistic at the end of Phase 2.
#'
#' \strong{Decision Rules:}
#'
#' * \strong{Phase 2 efficacy stop}: reject if the Phase-2 test statistic
#'   for the arm selected at rank \code{rankp0} satisfies
#'   \eqn{Z_{[rankp0]}(I_0) \ge b_0}.
#'
#' * \strong{Phase 2 futility stop}: stop for futility if the Phase-2 test
#'   statistic for the arm selected at rank \code{rankp0} satisfies
#'   \eqn{Z_{[rankp0]}(I_0) \le a_0}.
#'
#' * \strong{Continue to Phase 3}: if
#'   \eqn{a_0 < Z_{[rankp0]}(I_0) < b_0}, continue with the arm selected
#'   at rank \code{rankp0} only.
#'
#' * \strong{Phase 3 efficacy stop}: at look \eqn{k}, reject if the selected
#'   arm's Z-statistic exceeds the efficacy boundary and no earlier stop has
#'   occurred.
#'
#' * \strong{Phase 3 futility stop}: at look \eqn{k}, stop for futility if
#'   the selected arm's Z-statistic is below the futility boundary and no
#'   earlier stop has occurred.
#'
#' \strong{Design Assumptions:}
#'
#' * All active arms share the same information level in Phase 2.
#'
#' * Exactly one active arm is selected at the end of Phase 2 based on the
#'   \code{rankp0}-th largest observed Z-statistic when the trial continues
#'   to Phase 3.
#'
#' @return A list containing the following components:
#'
#' * \code{exitProbUpper}: A vector of length \eqn{K + 1}. The first element
#'   is the probability of stopping for efficacy in Phase 2; the remaining
#'   elements are the probabilities of stopping for efficacy at each look in
#'   Phase 3.
#'
#' * \code{exitProbLower}: A vector of length \eqn{K + 1}. The first element
#'   is the probability of stopping for futility in Phase 2; the remaining
#'   elements are the probabilities of stopping for futility at each look in
#'   Phase 3.
#'
#' * \code{exitProbByArmUpper}: A \eqn{(K + 1) \times M} matrix. The
#'   \eqn{(k, m)}-th entry gives the probability of stopping for efficacy at
#'   look \eqn{k} given that arm \eqn{m} is selected at rank \code{rankp0}.
#'
#' * \code{exitProbByArmLower}: A \eqn{(K + 1) \times M} matrix. The
#'   \eqn{(k, m)}-th entry gives the probability of stopping for futility at
#'   look \eqn{k} given that arm \eqn{m} is selected at rank \code{rankp0}.
#'
#' * \code{selectionProb}: A vector of length \eqn{M} containing the
#'   probability that each active arm is selected at rank \code{rankp0}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @examples
#'
#' # Setup: 2 active arms vs control in Phase 2; 1 selected arm vs control
#' # in Phase 3. Phase 3 has 2 sequential looks.
#'
#' # Information levels: equal spacing over 3 looks based on a maximum of
#' # 110 patients per arm, SD = 1.0
#' I <- c(110 / (2 * 1.0^2) * seq(1, 3)/3)
#'
#' # O'Brien-Fleming efficacy boundaries
#' b <- c(3.776605, 2.670463, 2.180424)
#'
#' # No futility stopping
#' p0 <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2,
#'                         b = b, I = I, nthreads = 1)
#' cumsum(p0$exitProbUpper)
#'
#' # Add futility stopping
#' a <- c(0, 0.5, b[3])
#' p1 <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2,
#'                         b = b, a = a, I = I, nthreads = 1)
#' cbind(
#'   cumulativeEfficacy = cumsum(p1$exitProbUpper),
#'   cumulativeFutility = cumsum(p1$exitProbLower)
#' )
#'
#' @export
exitprob_seamless <- function(M = NA_integer_,
                              r = 1,
                              theta = NA_real_,
                              corr_known = TRUE,
                              K = NA_integer_,
                              b = NULL,
                              a = NULL,
                              I = NULL,
                              rankp0 = 1L,
                              nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    exitprob_seamless_Rcpp(M, r, theta, corr_known, K, b, a, I, rankp0)
}


#' @title Efficacy Boundaries for Phase 2/3 Seamless Design
#' @description Calculates the efficacy stopping boundaries for a phase 2/3
#' seamless design, accounting for rank-based treatment selection
#' at the end of Phase 2 and sequential testing in Phase 3.
#'
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#' @param k The index of the current look in Phase 3.
#' @param informationRates A numeric vector of information rates up to the
#'   current look. Values must be strictly increasing and \eqn{\le 1}.
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A numeric vector of length \eqn{k+1} specifying the
#'   error spending time at each analysis. Values must be strictly increasing
#'   and \eqn{\le 1}. If omitted, defaults to \code{informationRates}.
#' @inheritParams param_efficacyStopping
#' @param rankp0 An integer between 1 and \code{M} specifying the rank of the
#'   arm to be selected at the end of Phase 2 for the purpose of determining
#'   the boundaries. For example, if \code{rankp0} is 1, the boundaries are
#'   determined based on the arm with the largest Z-statistic at the end of
#'   Phase 2; if \code{rankp0} is 2, the boundaries are determined based on the
#'   arm with the second largest Z-statistic at the end of Phase 2, and
#'   so on. The default is 1, which corresponds to the common practice of
#'   determining boundaries based on the top-ranked arm at the end of Phase 2.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' The function determines critical values by solving for the boundary that
#' satisfies the alpha-spending requirement, given the selection of the
#' arm at rank \code{rankp0} at the end of Phase 2.
#'
#' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
#' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
#' and \code{spendingTime} must be of full length \code{kMax}, and
#' \code{informationRates} and \code{spendingTime} must end with 1.
#'
#' @return A numeric vector of length \eqn{k + 1} containing the critical
#' values (on the standard normal Z-scale) for each analysis up to the
#' current look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @examples
#'
#' # Determine O'Brien-Fleming boundaries for a seamless design with
#' # 2 active arms in Phase 2 and 2 looks in Phase 3 (3 looks total).
#' getBound_seamless(M = 2, k = 2, informationRates = seq(1, 3)/3,
#'                   alpha = 0.025, typeAlphaSpending = "OF",
#'                   nthreads = 1)
#'
#' @export
getBound_seamless <- function(M = NA_integer_,
                              r = 1,
                              corr_known = TRUE,
                              k = NA_integer_,
                              informationRates = NA_real_,
                              alpha = 0.025,
                              typeAlphaSpending = "sfOF",
                              parameterAlphaSpending = NA_real_,
                              userAlphaSpending = NA_real_,
                              spendingTime = NA_real_,
                              efficacyStopping = NA_integer_,
                              rankp0 = 1L,
                              nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getBound_seamless_Rcpp(M, r, corr_known, k, informationRates,
                           alpha, typeAlphaSpending,
                           parameterAlphaSpending,
                           userAlphaSpending, spendingTime,
                           efficacyStopping, rankp0)
}


#' @title Power and Sample Size for Phase 2/3 Seamless Design
#' @description Computes either the maximum information and stopping
#' boundaries for a phase 2/3 seamless design, or the achieved power when
#' the maximum information and stopping boundaries are provided. Both
#' efficacy and futility stopping can be incorporated.
#'
#' @param beta Type II error rate. Provide either \code{beta} or \code{IMax};
#'   the other should be missing.
#' @param IMax Maximum information for any active arm versus the common
#'   control. Provide either \code{IMax} or \code{beta}; the other should
#'   be missing.
#' @param theta A vector of length \eqn{M} representing the true treatment
#'   effects for each active arm versus the common control. The global null
#'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
#'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#'   This option is only used for critical value calculations; the
#'   correlation is always derived from \eqn{r} for power calculations.
#' @param K Number of sequential looks in Phase 3.
#' @param informationRates A numeric vector of information rates fixed
#'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @param criticalValues The upper boundaries on the Z-statistic scale
#'   for the rank-selected arm in Phase 2 and the Z statistics for the
#'   selected arm in Phase 3.
#'   If missing, boundaries will be computed based on the specified alpha
#'   spending function.
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param futilityBounds A numeric vector of length \eqn{K} specifying
#'   futility boundaries on the Z scale at the end of Phase 2 for the
#'   rank-selected arm and on the Z scale for the \eqn{K - 1} analyses
#'   in Phase 3. The final analysis uses the efficacy boundary as the
#'   futility boundary.
#' @param futilityCP A numeric vector of length \eqn{K} specifying futility
#'   boundaries on the conditional power scale.
#' @param futilityTheta A numeric vector of length \eqn{K} specifying futility
#'   boundaries on the parameter scale.
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A numeric vector of length \eqn{K+1} specifying the
#'   error spending time at each analysis. Values must be strictly increasing
#'   and end at 1. If omitted, defaults to \code{informationRates}.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward. \code{rankp0 = 1} selects the largest
#'   Phase-2 Z-statistic, \code{rankp0 = 2} selects the second largest, and
#'   so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @return An S3 object of class \code{seamless} with the following components:
#'
#' * \code{overallResults}: A data frame containing:
#'     - \code{overallReject}: Overall probability of rejecting the null
#'       hypothesis.
#'     - \code{alpha}: Overall significance level.
#'     - \code{attainedAlpha}: The attained significance level, which may
#'       differ from \code{alpha} in the presence of futility stopping.
#'     - \code{M}: Number of active arms in Phase 2.
#'     - \code{r}: Randomization ratio per active arm versus control in
#'       Phase 2.
#'     - \code{corr_known}: Whether the phase-2 correlation was assumed known.
#'     - \code{rankp0}: The rank of the selected arm at the end of Phase 2.
#'     - \code{K}: Number of looks in Phase 3.
#'     - \code{information}: Maximum information for any active arm versus
#'       control.
#'     - \code{expectedInformationH1}: Expected information under the
#'       alternative.
#'     - \code{expectedInformationH0}: Expected information under the null.
#'     - \code{informationOverall}: Maximum information for the overall study.
#'     - \code{expectedInformationH1}: Expected information under the
#'       alternative for the overall study.
#'     - \code{expectedInformationH0}: Expected information under the null
#'       for the overall study.
#'
#' * \code{byStageResults}: A data frame containing:
#'     - \code{informationRates}: Information rates at each analysis.
#'     - \code{efficacyBounds}: Efficacy boundaries on the Z scale.
#'     - \code{futilityBounds}: Futility boundaries on the Z scale.
#'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
#'     - \code{futilityPerStage}: Probability of futility stopping at each
#'       stage.
#'     - \code{cumulativeRejection}: Cumulative probability of efficacy
#'       stopping.
#'     - \code{cumulativeFutility}: Cumulative probability of futility
#'       stopping.
#'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
#'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
#'     - \code{futilityTheta}: Futility boundaries on the parameter scale.
#'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
#'     - \code{futilityP}: Futility boundaries on the p-value scale.
#'     - \code{information}: Cumulative information at each analysis.
#'     - \code{informationOverall}: Cumulative information for the overall
#'       study at each analysis.
#'     - \code{efficacyStopping}: Indicator of whether efficacy stopping is
#'       permitted.
#'     - \code{futilityStopping}: Indicator of whether futility stopping is
#'       permitted.
#'     - \code{rejectPerStageH0}: Probability of efficacy stopping under the
#'       global null.
#'     - \code{futilityPerStageH0}: Probability of futility stopping under the
#'       global null.
#'     - \code{cumulativeRejectionH0}: Cumulative probability of efficacy
#'       stopping under the global null.
#'     - \code{cumulativeFutilityH0}: Cumulative probability of futility
#'       stopping under the global null.
#'
#' * \code{byArmResults}: A data frame containing:
#'     - \code{theta}: Parameter values for the active arms.
#'     - \code{selectionProb}: Probability an arm is selected at the
#'       end of Phase 2.
#'     - \code{powerByArm}: Probability of rejecting the null for each arm by
#'       trial end.
#'     - \code{condPowerByArm}: Conditional power for each arm given it was
#'       selected at rank \code{rankp0} at the end of Phase 2.
#'
#' * \code{settings}: A list of input settings:
#'     - \code{typeAlphaSpending}: Type of alpha spending function.
#'     - \code{parameterAlphaSpending}: Parameter value for the chosen alpha
#'       spending function.
#'     - \code{userAlphaSpending}: User-specified alpha spending values.
#'     - \code{typeBetaSpending}: Type of beta spending function.
#'     - \code{parameterBetaSpending}: Parameter value for the chosen beta
#'       spending function.
#'     - \code{userBetaSpending}: User-specified beta spending values.
#'     - \code{spendingTime}: Error-spending times at each analysis.
#'
#' @details
#' If \code{corr_known} is \code{FALSE}, critical boundaries are
#' computed assuming independence among the Phase-2 Wald statistics
#' (a conservative assumption when \code{rankp0 = 1}). Power calculations,
#' however, use the correlation implied by the randomization ratio \eqn{r}.
#'
#' Futility boundaries may be supplied directly on the Z scale, derived from
#' conditional power, derived from parameter values, or computed from a beta
#' spending function.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @examples
#'
#' # Example 1: obtain the maximum information given power with no futility
#' (design1 <- getDesign_seamless(
#'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
#'   K = 2, informationRates = seq(1, 3)/3,
#'   alpha = 0.025, typeAlphaSpending = "OF", nthreads = 1))
#'
#' # Example 2: obtain power given the maximum information and a futility rule
#' (design2 <- getDesign_seamless(
#'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
#'   K = 2, informationRates = seq(1, 3)/3,
#'   alpha = 0.025, typeAlphaSpending = "OF",
#'   futilityBounds = c(0.0, 0.5), nthreads = 1))
#'
#' # Example 3: derive futility boundaries using beta spending
#' (design3 <- getDesign_seamless(
#'   beta = 0.1, theta = c(-log(0.5), -log(0.7)),
#'   M = 2, r = 1.0, corr_known = FALSE,
#'   K = 2, informationRates = seq(1, 3)/3,
#'   alpha = 0.025, typeAlphaSpending = "sfOF",
#'   typeBetaSpending = "sfHSD", parameterBetaSpending = -2,
#'   nthreads = 1))
#'
#' @export
getDesign_seamless <- function(beta = NA_real_,
                               IMax = NA_real_,
                               theta = NA_real_,
                               M = NA_integer_,
                               r = 1,
                               corr_known = TRUE,
                               K = 1L,
                               informationRates = NA_real_,
                               efficacyStopping = NA_integer_,
                               futilityStopping = NA_integer_,
                               criticalValues = NULL,
                               alpha = 0.025,
                               typeAlphaSpending = "sfOF",
                               parameterAlphaSpending = NA_real_,
                               userAlphaSpending = NA_real_,
                               futilityBounds = NULL,
                               futilityCP = NULL,
                               futilityTheta = NULL,
                               typeBetaSpending = "none",
                               parameterBetaSpending = NA_real_,
                               userBetaSpending = NA_real_,
                               spendingTime = NA_real_,
                               rankp0 = 1L,
                               nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getDesign_seamless_Rcpp(beta, IMax, theta, M, r, corr_known,
                            K, informationRates, efficacyStopping,
                            futilityStopping, criticalValues,
                            alpha, typeAlphaSpending,
                            parameterAlphaSpending,
                            userAlphaSpending, futilityBounds,
                            futilityCP, futilityTheta,
                            typeBetaSpending,
                            parameterBetaSpending,
                            userBetaSpending, spendingTime,
                            rankp0)
}


#' @title Adaptive Phase 2/3 Seamless Design
#' @description
#' Calculates the conditional power for specified incremental
#' information, given the interim results, parameter value,
#' data-dependent changes in the error spending function, and
#' the number and spacing of interim looks. Conversely,
#' calculates the incremental information required to attain
#' a specified conditional power, given the interim results,
#' parameter value, data-dependent changes in
#' the error spending function, and the number and spacing of interim looks.
#'
#' @param betaNew The type II error for the secondary trial.
#' @param INew The maximum information for the active arm versus the common
#'   control in the secondary trial. Either
#'   \code{betaNew} or \code{INew} should be provided, while the other
#'   must be missing.
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#'   This option is only used for critical value calculations.
#' @param L The interim adaptation look in Phase 3.
#' @param zL The z-test statistic at the interim adaptation look of
#'   Phase 3.
#' @param theta The assumed treatment effect for the selected arm versus the
#'   common control.
#' @param IMax Maximum information for the active arm versus the common
#'   control for the original trial. Must be provided.
#' @param K Number of sequential looks in Phase 3.
#' @param informationRates A numeric vector of information rates fixed
#'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStopping Indicators of whether futility stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param criticalValues The upper boundaries on the z-test statistic
#'   scale for the rank-selected arm in Phase 2 and the z-test statistics
#'   for the selected arm in Phase 3 for the primary trial. If missing,
#'   boundaries will be computed based on the specified alpha spending function.
#' @param alpha The significance level of the primary trial.
#'   Defaults to 0.025.
#' @param typeAlphaSpending The type of alpha spending for the primary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
#'   \code{"user"} for user defined spending, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpending The parameter value of alpha spending
#'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param userAlphaSpending The user-defined alpha spending for the
#'   primary trial. Represents the cumulative alpha spent up to each stage.
#' @param futilityBounds The lower boundaries on the z-test statistic scale
#'   at the end of phase 2 for the rank-selected arm and on the z-test
#'   statistic scale in phase 3 for futility stopping for the primary trial.
#'   Defaults to \code{rep(-8, kMax-1)} if left unspecified.
#' @param futilityCP The conditional power-based futility bounds for the
#'   primary trial.
#' @param futilityTheta The parameter value-based futility bounds for the
#'   primary trial.
#' @param spendingTime The error spending time of the primary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRates}.
#' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
#'   for trial adaptation.
#' @param kNew The number of looks of the secondary trial.
#' @param informationRatesNew The spacing of looks of the secondary trial.
#' @param efficacyStoppingNew The indicators of whether efficacy stopping is
#'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStoppingNew The indicators of whether futility stopping is
#'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param typeAlphaSpendingNew The type of alpha spending for the secondary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param futilityBoundsInt The futility boundaries on the z statistic
#'   scale for new stages of the integrated trial.
#' @param futilityCPInt The conditional power-based futility bounds for
#'   new stages of the integrated trial.
#' @param futilityThetaInt The parameter value-based futility bounds for the
#'   new stages of the integrated trial.
#' @param typeBetaSpendingNew The type of beta spending for the secondary
#'   trial. One of the following:
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
#'   \code{"user"} for user defined spending, and
#'   \code{"none"} for no early futility stopping.
#'   Defaults to \code{"none"}.
#' @param parameterBetaSpendingNew The parameter value of beta spending
#'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
#'   and \eqn{\gamma} for \code{"sfHSD"}.
#' @param userBetaSpendingNew The user-defined cumulative beta spending.
#'   Represents the cumulative beta spent up to each stage of the
#'   secondary trial.
#' @param spendingTimeNew The error spending time of the secondary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRatesNew}.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward. \code{rankp0 = 1} selects the largest
#'   Phase-2 Z-statistic, \code{rankp0 = 2} selects the second largest, and
#'   so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @return An \code{adaptDesign_seamless} object with three list components:
#'
#' * \code{primaryTrial}: A list of selected information for the primary
#'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
#'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
#'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
#'   \code{information}, \code{alpha}, \code{conditionalAlpha},
#'   \code{conditionalPower}, and \code{MullerSchafer}.
#'
#' * \code{secondaryTrial}: A list of selected information for the seconary
#'   trial, including \code{overallReject}, \code{alpha}, \code{kMax},
#'   \code{maxInformation}, \code{informationRates}, \code{efficacyBounds},
#'   \code{futilityBounds}, \code{cumulativeRejection},
#'   \code{cumulativeFutility}, \code{cumulativeAlphaSpent},
#'   \code{information}, \code{typeAlphaSpending},
#'   \code{parameterAlphaSpending}, \code{typeBetaSpending},
#'   \code{parameterBetaSpending}, \code{userBetaSpending}, and
#'   \code{spendingTime}.
#'
#' * \code{integratedTrial}: A list of selected information for the integrated
#'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{K},
#'   \code{L}, \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
#'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
#'   and \code{information}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @seealso \code{\link{getDesign_seamless}}
#'
#' @examples
#'
#' (des1 <- adaptDesign_seamless(
#'   betaNew = 0.1, M = 2, r = 1, corr_known = FALSE,
#'   L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
#'   IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF", kNew = 1, nthreads = 1))
#'
#' @export
adaptDesign_seamless <- function(betaNew = NA_real_,
                                 INew = NA_real_,
                                 M = NA_integer_,
                                 r = 1,
                                 corr_known = TRUE,
                                 L = NA_integer_,
                                 zL = NA_real_,
                                 theta = NA_real_,
                                 IMax = NA_real_,
                                 K = NA_integer_,
                                 informationRates = NA_real_,
                                 efficacyStopping = NA_integer_,
                                 futilityStopping = NA_integer_,
                                 criticalValues = NULL,
                                 alpha = 0.025,
                                 typeAlphaSpending = "sfOF",
                                 parameterAlphaSpending = NA_real_,
                                 userAlphaSpending = NA_real_,
                                 futilityBounds = NULL,
                                 futilityCP = NULL,
                                 futilityTheta = NULL,
                                 spendingTime = NA_real_,
                                 MullerSchafer = FALSE,
                                 kNew = NA_integer_,
                                 informationRatesNew = NA_real_,
                                 efficacyStoppingNew = NA_integer_,
                                 futilityStoppingNew = NA_integer_,
                                 typeAlphaSpendingNew = "sfOF",
                                 parameterAlphaSpendingNew = NA_real_,
                                 futilityBoundsInt = NULL,
                                 futilityCPInt = NULL,
                                 futilityThetaInt = NULL,
                                 typeBetaSpendingNew = "none",
                                 parameterBetaSpendingNew = NA_real_,
                                 userBetaSpendingNew = NA_real_,
                                 spendingTimeNew = NA_real_,
                                 rankp0 = 1L,
                                 nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    adaptDesign_seamless_Rcpp(betaNew, INew, M, r, corr_known,
                              L, zL, theta, IMax, K, informationRates,
                              efficacyStopping, futilityStopping,
                              criticalValues, alpha, typeAlphaSpending,
                              parameterAlphaSpending,
                              userAlphaSpending, futilityBounds,
                              futilityCP, futilityTheta,
                              spendingTime, MullerSchafer,
                              kNew, informationRatesNew,
                              efficacyStoppingNew,
                              futilityStoppingNew,
                              typeAlphaSpendingNew,
                              parameterAlphaSpendingNew,
                              futilityBoundsInt, futilityCPInt,
                              futilityThetaInt, typeBetaSpendingNew,
                              parameterBetaSpendingNew,
                              userBetaSpendingNew, spendingTimeNew,
                              rankp0)
}


#' @title Conditional Power for Phase 2/3 Seamless Design
#' @description Obtains the conditional power for specified incremental
#' information given the interim results, parameter values, and
#' data-dependent changes in the error spending function, as well as the
#' number and spacing of interim looks.
#'
#' @param INew The maximum information for the active arm versus the common
#'   control in the secondary trial.
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#'   This option is only used for critical value calculations.
#' @param L The interim adaptation look in Phase 3.
#' @param zL The z-test statistic at the interim adaptation look of
#'   Phase 3.
#' @param theta The assumed treatment effect for the selected arm versus the
#'   common control.
#' @param IMax Maximum information for the active arm versus the common
#'   control for the original trial. Must be provided.
#' @param K Number of sequential looks in Phase 3.
#' @param informationRates A numeric vector of information rates fixed
#'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStopping Indicators of whether futility stopping is
#'   allowed at each stage of the primary trial. Defaults to true
#'   if left unspecified.
#' @param criticalValues The upper boundaries on the z-test statistic scale
#'   for the rank-selected arm in Phase 2 and the z-test statistics for the
#'   selected arm in Phase 3 for the primary trial. If missing, boundaries
#'   will be computed based on the specified alpha spending function.
#' @param alpha The significance level of the primary trial.
#'   Defaults to 0.025.
#' @param typeAlphaSpending The type of alpha spending for the primary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
#'   \code{"user"} for user defined spending, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpending The parameter value of alpha spending
#'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param userAlphaSpending The user-defined alpha spending for the
#'   primary trial. Represents the cumulative alpha spent up to each stage.
#' @param futilityBounds The lower boundaries on the z-test statistic
#'   scale for the rank-selected arm in Phase 2 and the z-test statistics
#'   for the selected arm in Phase 3 for the primary trial.
#' @param futilityCP The conditional power-based futility bounds for the
#'   primary trial.
#' @param futilityTheta The parameter value-based futility bounds for the
#'   primary trial.
#' @param spendingTime The error spending time of the primary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRates}.
#' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
#'   for trial adaptation.
#' @param kNew The number of looks of the secondary trial.
#' @param informationRatesNew The spacing of looks of the secondary trial.
#' @param efficacyStoppingNew The indicators of whether efficacy stopping is
#'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStoppingNew The indicators of whether futility stopping is
#'   allowed at each look of the secondary trial. Defaults to true
#'   if left unspecified.
#' @param typeAlphaSpendingNew The type of alpha spending for the secondary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param futilityBoundsInt The futility boundaries on the z statistic
#'   scale for new stages of the integrated trial.
#' @param futilityCPInt The conditional power-based futility bounds for
#'   new stages of the integrated trial.
#' @param futilityThetaInt The parameter value-based futility bounds for the
#'   new stages of the integrated trial.
#' @param typeBetaSpendingNew The type of beta spending for the secondary
#'   trial. One of the following:
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early futility stopping.
#'   Defaults to \code{"none"}.
#' @param parameterBetaSpendingNew The parameter value of beta spending
#'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
#'   and \eqn{\gamma} for \code{"sfHSD"}.
#' @param spendingTimeNew The error spending time of the secondary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRatesNew}.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward when the trial continues to Phase 3.
#'   \code{rankp0 = 1} selects the largest Phase-2 Z-statistic,
#'   \code{rankp0 = 2} selects the second largest, and so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @return A vector of two conditional powers given the interim results and
#' parameter values, one without design change and the other with
#' data-dependent design changes.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @seealso \code{\link{adaptDesign_seamless}}
#'
#' @examples
#'
#' getCP_seamless(
#'   INew = 198 / 4, M = 2, r = 1, corr_known = FALSE,
#'   L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
#'   IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF", kNew = 1,
#'   nthreads = 1)
#'
#' @export
getCP_seamless <- function(INew = NA_real_,
                           M = NA_integer_,
                           r = 1,
                           corr_known = TRUE,
                           L = NA_integer_,
                           zL = NA_real_,
                           theta = NA_real_,
                           IMax = NA_real_,
                           K = NA_integer_,
                           informationRates = NA_real_,
                           efficacyStopping = NA_integer_,
                           futilityStopping = NA_integer_,
                           criticalValues = NULL,
                           alpha = 0.025,
                           typeAlphaSpending = "sfOF",
                           parameterAlphaSpending = NA_real_,
                           userAlphaSpending = NA_real_,
                           futilityBounds = NULL,
                           futilityCP = NULL,
                           futilityTheta = NULL,
                           spendingTime = NA_real_,
                           MullerSchafer = FALSE,
                           kNew = NA_integer_,
                           informationRatesNew = NA_real_,
                           efficacyStoppingNew = NA_integer_,
                           futilityStoppingNew = NA_integer_,
                           typeAlphaSpendingNew = "sfOF",
                           parameterAlphaSpendingNew = NA_real_,
                           futilityBoundsInt = NULL,
                           futilityCPInt = NULL,
                           futilityThetaInt = NULL,
                           typeBetaSpendingNew = "none",
                           parameterBetaSpendingNew = NA_real_,
                           spendingTimeNew = NA_real_,
                           rankp0 = 1L,
                           nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getCP_seamless_Rcpp(INew, M, r, corr_known, L, zL, theta, IMax,
                        K, informationRates, efficacyStopping,
                        futilityStopping, criticalValues,
                        alpha, typeAlphaSpending,
                        parameterAlphaSpending, userAlphaSpending,
                        futilityBounds, futilityCP, futilityTheta,
                        spendingTime, MullerSchafer,
                        kNew, informationRatesNew,
                        efficacyStoppingNew, futilityStoppingNew,
                        typeAlphaSpendingNew,
                        parameterAlphaSpendingNew,
                        futilityBoundsInt, futilityCPInt,
                        futilityThetaInt, typeBetaSpendingNew,
                        parameterBetaSpendingNew, spendingTimeNew,
                        rankp0)
}


#' @title Confidence Interval After Trial Termination for Phase 2/3
#' Seamless Design
#' @description Obtains the p-value, point estimate, and
#' confidence interval after the end of a phase 2/3 seamless trial.
#'
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#'   This option is only used for critical value calculations.
#' @param L The termination look in Phase 3.
#' @param zL The z-test statistic at the termination look.
#' @param IMax Maximum information for any active arm versus the common
#'   control.
#' @param informationRates The information rates up to look \code{L}.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage up to look \code{L}.
#'   Defaults to \code{TRUE} if left unspecified.
#' @param criticalValues The upper boundaries on the z-test statistic
#'   scale for the rank-selected arm in Phase 2 and the z-test statistics
#'   for the selected arm
#'   in Phase 3 up to look \code{L}. If missing, boundaries will be
#'   computed based on the specified alpha spending function.
#' @inheritParams param_alpha
#' @param typeAlphaSpending The type of alpha spending for the trial.
#'   One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @inheritParams param_parameterAlphaSpending
#' @param spendingTime The error spending time up to look \code{L}.
#'   Defaults to missing, in which case, it is the same as
#'   \code{informationRates}.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward. \code{rankp0 = 1} selects the largest
#'   Phase-2 Z-statistic, \code{rankp0 = 2} selects the second largest, and
#'   so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
#' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
#' and \code{spendingTime} must be of full length \eqn{K + 1}, and
#' \code{informationRates} and \code{spendingTime} must end with 1.
#'
#' @return A data frame with the following components:
#'
#' * \code{pvalue}: p-value for rejecting the null hypothesis.
#'
#' * \code{thetahat}: Point estimate of the parameter.
#'
#' * \code{cilevel}: Confidence interval level.
#'
#' * \code{lower}: Lower bound of confidence interval.
#'
#' * \code{upper}: Upper bound of confidence interval.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive two-stage seamless sequential design for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
#'
#' @examples
#' getCI_seamless(
#'   L = 2, zL = 2.075, M = 2, r = 1, corr_known = FALSE,
#'   IMax = 300 / 4, informationRates = c(1/3, 2/3, 1),
#'   alpha = 0.025, typeAlphaSpending = "sfOF", nthreads = 1)
#'
#' @export
getCI_seamless <- function(M = NA_integer_,
                           r = 1,
                           corr_known = TRUE,
                           L = NA_integer_,
                           zL = NA_real_,
                           IMax = NA_real_,
                           informationRates = NA_real_,
                           efficacyStopping = NA_integer_,
                           criticalValues = NA_real_,
                           alpha = 0.025,
                           typeAlphaSpending = "sfOF",
                           parameterAlphaSpending = NA_real_,
                           spendingTime = NA_real_,
                           rankp0 = 1L,
                           nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
  getCI_seamless_Rcpp(M, r, corr_known, L, zL, IMax, informationRates,
                      efficacyStopping, criticalValues, alpha,
                      typeAlphaSpending, parameterAlphaSpending,
                      spendingTime, rankp0)
}


#' @title Confidence Interval After Adaptation for Phase 2/3 Seamless
#' Design
#' @description Obtains the p-value, conservative point estimate, and
#' confidence interval after the end of an adaptive phase 2/3 seamless
#' design.
#'
#' @param M Number of active treatment arms in Phase 2.
#' @param r Randomization ratio of each active arm to the common control
#'   in Phase 2.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics in Phase 2 is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used, which is only valid when \code{rankp0 = 1} (i.e., the arm
#'   with the largest Phase-2 Z-statistic is selected for Phase 3).
#'   This option is only used for critical value calculations.
#' @param L The interim adaptation look in Phase 3.
#' @param zL The z-test statistic at the interim adaptation look of
#'   Phase 3.
#' @param IMax Maximum information for the active arm versus the common
#'   control for the original trial. Must be provided.
#' @param K Number of sequential looks in Phase 3.
#' @param informationRates A numeric vector of information rates fixed
#'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param criticalValues The upper boundaries on the z-test statistic
#'   scale for the rank-selected arm in Phase 2 and the z-test statistics
#'   for the selected arm
#'   in Phase 3 for the primary trial. If missing, boundaries
#'   will be computed based on the specified alpha spending function.
#' @param alpha The significance level of the primary trial.
#'   Defaults to 0.025.
#' @param typeAlphaSpending The type of alpha spending for the primary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpending The parameter value of alpha spending
#'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param spendingTime The error spending time of the primary trial.
#'   Defaults to missing, in which case, it is the same as
#'   \code{informationRates}.
#' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
#'   for trial adaptation.
#' @param Lc The termination look of the integrated trial.
#' @param zLc The z-test statistic at the termination look of the
#'   integrated trial.
#' @param INew The maximum information for the active arm versus the common
#'   control in the secondary trial.
#' @param informationRatesNew The spacing of looks of the secondary trial.
#' @param efficacyStoppingNew The indicators of whether efficacy stopping is
#'   allowed at each look of the secondary trial.
#'   Defaults to \code{TRUE} if left unspecified.
#' @param typeAlphaSpendingNew The type of alpha spending for the secondary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"P"} for Pocock boundaries,
#'   \code{"WT"} for Wang & Tsiatis boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param spendingTimeNew The error spending time of the secondary trial.
#'   Defaults to missing, in which case, it is
#'   the same as \code{informationRatesNew}.
#' @param rankp0 An integer between 1 and \code{M} specifying which ranked
#'   Phase-2 arm is carried forward. \code{rankp0 = 1} selects the largest
#'   Phase-2 Z-statistic, \code{rankp0 = 2} selects the second largest, and
#'   so on.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' If typeAlphaSpendingNew is \code{"OF"}, \code{"P"}, \code{"WT"},
#' or \code{"none"}, then \code{informationRatesNew},
#' \code{efficacyStoppingNew}, and \code{spendingTimeNew} must be of full
#' length \code{kNew}, and \code{informationRatesNew} and
#' \code{spendingTimeNew} must end with 1.
#'
#' @return A data frame with the following variables:
#'
#' * \code{pvalue}: p-value for rejecting the null hypothesis.
#'
#' * \code{thetahat}: Point estimate of the parameter.
#'
#' * \code{cilevel}: Confidence interval level.
#'
#' * \code{lower}: Lower bound of confidence interval.
#'
#' * \code{upper}: Upper bound of confidence interval.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @examples
#' getADCI_seamless(
#'   M = 2, r = 1, corr_known = FALSE,
#'   L = 1, zL = -log(0.67) * sqrt(80 / 4),
#'   IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF",
#'   Lc = 2, zLc = -log(0.677) * sqrt(236 / 4), INew = 236 / 4,
#'   nthreads = 1)
#'
#' @export
getADCI_seamless <- function(M = NA_integer_,
                             r = 1,
                             corr_known = TRUE,
                             L = NA_integer_,
                             zL = NA_real_,
                             IMax = NA_real_,
                             K = NA_integer_,
                             informationRates = NA_real_,
                             efficacyStopping = NA_integer_,
                             criticalValues = NA_real_,
                             alpha = 0.25,
                             typeAlphaSpending = "sfOF",
                             parameterAlphaSpending = NA_real_,
                             spendingTime = NA_real_,
                             MullerSchafer = FALSE,
                             Lc = NA_integer_,
                             zLc = NA_real_,
                             INew = NA_real_,
                             informationRatesNew = NA_real_,
                             efficacyStoppingNew = NA_integer_,
                             typeAlphaSpendingNew = "sfOF",
                             parameterAlphaSpendingNew = NA_real_,
                             spendingTimeNew = NA_real_,
                             rankp0 = 1L,
                             nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
  getADCI_seamless_Rcpp(M, r, corr_known, L, zL, IMax, K,
                        informationRates, efficacyStopping,
                        criticalValues, alpha, typeAlphaSpending,
                        parameterAlphaSpending, spendingTime,
                        MullerSchafer, Lc, zLc, INew,
                        informationRatesNew, efficacyStoppingNew,
                        typeAlphaSpendingNew,
                        parameterAlphaSpendingNew,
                        spendingTimeNew, rankp0)
}
