#' @title Exit Probabilities for Multi-Arm Multi-Stage Design
#' @description Computes the exit (rejection) probabilities for a multi-arm
#' multi-stage design.
#'
#' @param M Number of active treatment arms.
#' @param r Randomization ratio of each active arm to the common control.
#' @param theta A vector of length \eqn{M} representing the true treatment
#'   effects for each active arm versus the common control.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used.
#' @param kMax Number of sequential looks.
#' @param b A vector of efficacy boundaries for the max-Z statistics.
#' @param a A vector of futility boundaries for the max-Z statistics.
#' @param I A vector of information levels for any active arm versus the
#'   common control.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' The function assumes a multivariate normal distribution for the Wald
#' statistics and all active arms share the same information level.
#'
#' @return A vector \code{exitProb} of length \code{kMax} containing the
#' probability of rejection at each look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @examples
#'
#' # Setup: 2 active arms vs control and 3 sequential looks.
#'
#' # Information levels: equal spacing over 3 looks based on a maximum of
#' # 95 patients per arm, SD = 1.0
#' I <- 95 / (2 * 1.0^2) * seq(1, 3)/3
#'
#' # O'Brien-Fleming critical values
#' b <- c(3.886562, 2.748214, 2.243907)
#'
#' # Type I error under the global null hypothesis
#' p0 <- exitprob_multiarm(M = 2, theta = c(0, 0), kMax = 3,
#'                         b = b, I = I, nthreads = 1)
#' cumsum(p0$exitProbUpper)
#'
#' # Power under alternative: Treatment effects of 0.3 and 0.5
#' p1 <- exitprob_multiarm(M = 2, theta = c(0.3, 0.5), kMax = 3,
#'                         b = b, I = I, nthreads = 1)
#' cumsum(p1$exitProbUpper)
#'
#' @export
exitprob_multiarm <- function(M = NA_integer_,
                              r = 1,
                              theta = NA_real_,
                              corr_known = TRUE,
                              kMax = NA_integer_,
                              b = NULL,
                              a = NULL,
                              I = NULL,
                              nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    exitprob_multiarm_Rcpp(M, r, theta, corr_known, kMax, b, a, I)
}


#' @title Efficacy Boundaries for Multi-Arm Multi-Stage Design
#' @description Calculates the efficacy stopping boundaries for a multi-arm
#' multi-stage design.
#'
#' @param M Number of active treatment arms.
#' @param r Randomization ratio of each active arm to the common control.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is assumed.
#' @param k The index of the current look.
#' @param informationRates A numeric vector of information rates up to the
#'   current look. Values must be strictly increasing and \eqn{\le 1}.
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param spendingTime A numeric vector of length \eqn{k} specifying the
#'   error spending time at each analysis. Values must be strictly increasing
#'   and \eqn{\le 1}. If omitted, defaults to \code{informationRates}.
#' @inheritParams param_efficacyStopping
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' The function determines critical values by solving for the boundary that
#' satisfies the alpha-spending requirement.
#'
#' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
#' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
#' and \code{spendingTime} must be of full length \code{kMax}, and
#' \code{informationRates} and \code{spendingTime} must end with 1.
#'
#' @return A numeric vector of length \eqn{k} containing the critical
#' values (on the standard normal Z-scale) for each analysis up to the
#' current look.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @examples
#'
#' # Determine O'Brien-Fleming boundaries for a TSSSD with
#' # 2 active arms and 3 looks.
#' getBound_multiarm(M = 2, k = 3, informationRates = seq(1, 3)/3,
#'                   alpha = 0.025, typeAlphaSpending = "OF",
#'                   nthreads = 1)
#'
#' @export
getBound_multiarm <- function(M = NA_integer_,
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
                              nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getBound_multiarm_Rcpp(M, r, corr_known, k, informationRates, alpha,
                           typeAlphaSpending, parameterAlphaSpending,
                           userAlphaSpending, spendingTime,
                           efficacyStopping)
}


#' @title Power and Sample Size for Multi-Arm Multi-Stage Design
#' @description Computes either the maximum information and stopping
#' boundaries for a multi-arm multi-stage design, or
#' the achieved power when the maximum information and stopping boundaries
#' are provided.
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
#' @param M Number of active treatment arms.
#' @param r Randomization ratio of each active arm to the common control.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used.
#' @param kMax Number of sequential looks.
#' @param informationRates A numeric vector of information rates fixed
#'   before the trial. If unspecified, defaults to \eqn{(1:kMax) / kMax}.
#' @inheritParams param_efficacyStopping
#' @inheritParams param_futilityStopping
#' @param criticalValues The matrix of by-level upper boundaries on the
#'   max z-test statistic scale for efficacy stopping.
#'   The first column is for level \code{M}, the second column is for
#'   level \code{M - 1}, and so on, with the last column for level 1.
#'   If left unspecified, the critical values will be computed based
#'   on the specified alpha spending function.
#' @inheritParams param_alpha
#' @inheritParams param_typeAlphaSpending
#' @inheritParams param_parameterAlphaSpending
#' @inheritParams param_userAlphaSpending
#' @param futilityBounds A numeric vector of length \code{kMax - 1}
#'   specifying the futility boundaries on the max z-test statistic scale
#'   for futility stopping.
#' @param futilityCP A numeric vector of length \code{kMax - 1} specifying
#'   the futility boundaries on the conditional power scale for futility
#'   stopping.
#' @param futilityTheta A numeric vector of length \code{kMax - 1} specifying
#'   the futility boundaries on the parameter scale for futility stopping.
#' @inheritParams param_typeBetaSpending
#' @inheritParams param_parameterBetaSpending
#' @inheritParams param_userBetaSpending
#' @param spendingTime A numeric vector of length \code{kMax} specifying the
#'   error spending time at each analysis. Values must be strictly increasing
#'   and ends at 1. If omitted, defaults to \code{informationRates}.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @return An S3 object of class \code{multiarm} with the following components:
#'
#' * \code{overallResults}: A data frame containing:
#'     - \code{overallReject}: Overall probability of rejecting the global
#'       null hypothesis.
#'     - \code{alpha}: Overall significance level.
#'     - \code{attainedAlpha}: The attained significance level, which is
#'       different from the overall significance level in the presence of
#'       futility stopping.
#'     - \code{M}: Number of active arms.
#'     - \code{r}: Randomization ratio per active arm versus control.
#'     - \code{corr_known}: Whether the correlation among Wald statistics
#'       was assumed known.
#'     - \code{kMax}: Number of stages.
#'     - \code{information}: Maximum information for any active arm versus
#'       control.
#'     - \code{expectedInformationH1}: The expected information under H1.
#'     - \code{expectedInformationH0}: The expected information under H0.
#'
#' * \code{byStageResults}: A data frame containing:
#'     - \code{informationRates}: Information rates at each analysis.
#'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale.
#'     - \code{futilityBounds}: Futility boundaries on the max Z-scale.
#'     - \code{rejectPerStage}: Probability of efficacy stopping at each stage.
#'     - \code{futilityPerStage}: Probability of futility stopping at each stage.
#'     - \code{cumulativeRejection}: Cumulative probability of efficacy stopping.
#'     - \code{cumulativeFutility}: Cumulative probability of futility stopping.
#'     - \code{cumulativeAlphaSpent}: Cumulative alpha spent.
#'     - \code{efficacyTheta}: Efficacy boundaries on the parameter scale.
#'     - \code{futilityTheta}: Futility boundaries on the parameter scale.
#'     - \code{efficacyP}: Efficacy boundaries on the p-value scale.
#'     - \code{futilityP}: Futility boundaries on the p-value scale.
#'     - \code{information}: Cumulative information for any active arm versus
#'       control at each analysis.
#'     - \code{efficacyStopping}: Indicator of whether efficacy stopping
#'       is permitted at each stage.
#'     - \code{futilityStopping}: Indicator of whether futility stopping
#'       is permitted at each stage.
#'     - \code{rejectPerStageH0}: The probability for efficacy stopping
#'       under H0.
#'     - \code{futilityPerStageH0}: The probability for futility stopping
#'       under H0.
#'     - \code{cumulativeRejectionH0}: The cumulative probability for
#'       efficacy stopping under H0.
#'     - \code{cumulativeFutilityH0}: The cumulative probability for
#'       futility stopping under H0.
#'
#' * \code{settings}: A list of input settings:
#'     - \code{typeAlphaSpending}: The type of alpha spending.
#'     - \code{parameterAlphaSpending}: The parameter value for the chosen
#'       alpha spending function.
#'     - \code{userAlphaSpending}: The user-specified alpha spending values.
#'     - \code{typeBetaSpending}: The type of beta spending.
#'     - \code{parameterBetaSpending}: The parameter value for the chosen
#'       beta spending function.
#'     - \code{userBetaSpending}: The user-specified beta spending values.
#'     - \code{spendingTime}: The error-spending time at each analysis.
#'
#' * \code{byLevelBounds}: A data frame containing the efficacy boundaries
#'   for each level of testing (i.e., number of active arms remaining) and
#'   each stage. Columns include:
#'     - \code{level}: Number of active arms remaining (1 to \eqn{M}).
#'     - \code{stage}: Stage index (1 to \code{kMax}).
#'     - \code{efficacyBounds}: Efficacy boundaries on the max Z-scale
#'       for the given level and stage.
#'
#' @details If \code{corr_known} is \code{FALSE}, critical boundaries are
#' computed assuming independence among the Wald statistics in each stage
#' (a conservative assumption). Power calculations, however, use the
#' correlation implied by the randomization ratio \eqn{r}.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @examples
#'
#' # Example 1: obtain the maximum information given power
#' (design1 <- getDesign_multiarm(
#'   beta = 0.1, theta = c(0.3, 0.5), M = 2, r = 1.0,
#'   kMax = 3, informationRates = seq(1, 3)/3,
#'   alpha = 0.025, typeAlphaSpending = "OF", nthreads = 1))
#'
#' # Example 2: obtain power given the maximum information
#' (design2 <- getDesign_multiarm(
#'   IMax = 110/(2*1^2), theta = c(0.3, 0.5), M = 2, r = 1.0,
#'   kMax = 3, informationRates = seq(1, 3)/3,
#'   alpha = 0.025, typeAlphaSpending = "OF", nthreads = 1))
#'
#' @export
getDesign_multiarm <- function(beta = NA_real_,
                               IMax = NA_real_,
                               theta = NA_real_,
                               M = NA_integer_,
                               r = 1,
                               corr_known = TRUE,
                               kMax = 1L,
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
                               nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getDesign_multiarm_Rcpp(beta, IMax, theta, M, r, corr_known,
                            kMax, informationRates, efficacyStopping,
                            futilityStopping, criticalValues,
                            alpha, typeAlphaSpending,
                            parameterAlphaSpending, userAlphaSpending,
                            futilityBounds, futilityCP,
                            futilityTheta, typeBetaSpending,
                            parameterBetaSpending,
                            userBetaSpending, spendingTime)
}


#' @title Adaptive Multi-Arm Multi-Stage Design
#' @description
#' Calculates the conditional power for specified incremental
#' information, given the interim results, parameter value,
#' data-dependent changes in treatment selection,
#' the error spending function, and
#' the number and spacing of interim looks. Conversely,
#' calculates the incremental information required to attain
#' a specified conditional power, given the interim results,
#' parameter value, data-dependent changes in treatment selection,
#' the error spending function, and the number and spacing of interim looks.
#'
#' @param betaNew The type II error for the secondary trial.
#' @param INew The maximum information for any active arm versus the common
#'   control in the secondary trial. Either
#'   \code{betaNew} or \code{INew} should be provided, while the other
#'   must be missing.
#' @param M Number of active treatment arms in the primary trial.
#' @param r Randomization ratio of each active arm to the common control
#'   in the primary trial.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is assumed.
#' @param L The interim adaptation look of the primary trial.
#' @param zL The z-test statistics at the interim adaptation look of
#'   the primary trial.
#' @param theta A vector of length \eqn{M} representing the assumed treatment
#'   effects for each active arm versus the common control. The global null
#'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
#'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
#' @param IMax Maximum information for any active arm versus the common
#'   control for the primary trial. Must be provided.
#' @param kMax The maximum number of stages of the primary trial.
#' @param informationRates The information rates of the primary trial.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStopping Indicators of whether futility stopping is
#'   allowed at each stage of the primary trial.
#' @param criticalValues The matrix of by-level upper boundaries on the
#'   max z-test statistic scale for efficacy stopping for the primary trial.
#'   The first column is for level \code{M}, the second column is for
#'   level \code{M - 1}, and so on, with the last column for level 1.
#'   If left unspecified, the critical values will be computed based
#'   on the specified alpha spending function.
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
#' @param futilityBounds The futility boundaries on the max-z statistic
#'   scale for the primary trial. Defaults to \code{rep(-8, kMax-1)}
#'   if left unspecified.
#' @param futilityCP The conditional power-based futility bounds for the
#'   primary trial.
#' @param futilityTheta The parameter value-based futility bounds for the
#'   primary trial.
#' @param spendingTime The error spending time of the primary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRates}.
#' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
#'   for trial adaptation.
#' @param MNew Number of active treatment arms in the secondary trial.
#' @param selected The indices of the selected active treatment arms for
#'   the secondary trial.
#' @param rNew Randomization ratio of each active arm to the common control
#'   in the secondary trial.
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
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param futilityBoundsInt The futility boundaries on the max-z statistic
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
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @return An \code{adaptDesign_multiarm} object with three list components:
#'
#' * \code{primaryTrial}: A list of selected information for the primary
#'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{L},
#'   \code{zL}, \code{theta}, \code{maxInformation}, \code{kMax},
#'   \code{informationRates}, \code{efficacyBounds}, \code{futilityBounds},
#'   \code{information}, \code{alpha}, \code{conditionalAlpha},
#'   \code{conditionalPower}, \code{MullerSchafer}, and \code{byLevelBounds},
#'   where \code{byLevelBounds} is a data frame with columns \code{level},
#'   \code{stage}, and \code{efficacyBounds}, representing the efficacy
#'   bounds for each combination of the number of active arms and
#'   the stage of analysis in the primary trial.
#'
#' * \code{secondaryTrial}: A list of selected information for the secondary
#'   trial, including \code{overallReject}, \code{alpha}, \code{M}, \code{r},
#'   \code{selected}, \code{corr_known}, \code{kMax}, \code{maxInformation},
#'   \code{informationRates}, \code{cumulativeRejection},
#'   \code{cumulativeAlphaSpent}, \code{information},
#'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
#'   \code{typeBetaSpending}, \code{parameterBetaSpending},
#'   \code{userBetaSpending}, \code{spendingTime}, and
#'   \code{byHypothesisBounds}, where \code{byHypothesisBounds} is a
#'   data frame with columns \code{hypothesis}, \code{stage},
#'   \code{efficacyBounds}, and \code{futilityBounds}, representing
#'   the efficacy and futility bounds for each hypothesis and each
#'   stage of analysis in the secondary trial.
#'
#' * \code{integratedTrial}: A list of selected information for the integrated
#'   trial, including \code{M}, \code{r}, \code{corr_known}, \code{MNew},
#'   \code{rNew}, \code{selected}, \code{L}, \code{zL}, \code{theta},
#'   \code{maxInformation}, \code{kMax}, \code{informationRates},
#'   \code{efficacyBounds}, \code{futilityBounds}, \code{information},
#'   and \code{byIntersectionBounds}, where \code{byIntersectionBounds} is
#'   a data frame with columns \code{intersectionHypothesis}, \code{stage},
#'   and \code{efficacyBounds}, representing the efficacy bounds for
#'   each intersection hypothesis and each stage of analysis in the
#'   integrated trial.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Ping Gao, Yingqiu Li.
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @seealso \code{\link{getDesign_multiarm}}
#'
#' @examples
#'
#' # Two active treatment arms are compared with a common control in a
#' # two-look time-to-event design using O'Brien–Fleming–type alpha spending.
#' # Suppose each active arm has a true hazard ratio of 0.75 versus control,
#' # and the total number of events across all three arms at the final analysis
#' # is 486. This corresponds to approximately 324 events for each active arm
#' # versus the common control. Under these assumptions, the trial has about
#' # 80% power to detect the treatment effect in at least one active arm.
#'
#' (des1 <- getDesign_multiarm(
#'   IMax = 324 / 4, theta = c(-log(0.75), -log(0.75)),
#'   M = 2, r = 1, kMax = 2, informationRates = c(1/2, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF", nthreads = 1))
#'
#' # Now assume that, at the interim analysis, the observed hazard ratios for
#' # the two active arms versus control are 0.91 and 0.78, respectively. Using
#' # the rule “drop any arm with an observed hazard ratio > 0.9”, arm 1 is
#' # dropped. We then aim to achieve 80% conditional power to detect a hazard
#' # ratio of 0.78 for the remaining arm at the final look. The analysis below
#' # indicates that the required total number of events for arm 2 versus control
#' # at the final analysis should be increased from 324 to 535.
#'
#' (des2 <- adaptDesign_multiarm(
#'   betaNew = 0.2, M = 2, r = 1, corr_known = FALSE,
#'   L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
#'   theta = c(-log(0.91), -log(0.78)),
#'   IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF",
#'   MNew = 1, selected = 2, rNew = 1, nthreads = 1))
#'
#' @export
adaptDesign_multiarm <- function(betaNew = NA_real_,
                                 INew = NA_real_,
                                 M = NA_integer_,
                                 r = 1,
                                 corr_known = TRUE,
                                 L = NA_integer_,
                                 zL = NA_real_,
                                 theta = NA_real_,
                                 IMax = NA_real_,
                                 kMax = NA_integer_,
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
                                 MNew = NA_integer_,
                                 selected = NA_integer_,
                                 rNew = 1,
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
                                 nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    adaptDesign_multiarm_Rcpp(betaNew, INew, M, r, corr_known,
                              L, zL, theta, IMax, kMax,
                              informationRates, efficacyStopping,
                              futilityStopping, criticalValues,
                              alpha, typeAlphaSpending,
                              parameterAlphaSpending,
                              userAlphaSpending,
                              futilityBounds, futilityCP,
                              futilityTheta, spendingTime,
                              MullerSchafer, MNew, selected,
                              rNew, kNew, informationRatesNew,
                              efficacyStoppingNew,
                              futilityStoppingNew,
                              typeAlphaSpendingNew,
                              parameterAlphaSpendingNew,
                              futilityBoundsInt, futilityCPInt,
                              futilityThetaInt, typeBetaSpendingNew,
                              parameterBetaSpendingNew,
                              userBetaSpendingNew, spendingTimeNew)
}


#' @title Conditional Power for Multi-Arm Multi-Stage Design
#' @description Obtains the conditional power for specified incremental
#' information given the interim results, parameter values, and
#' data-dependent changes in the selected treatment(s),
#' the error spending function, as well as the
#' number and spacing of interim looks.
#'
#' @param INew The maximum information for any active arm versus the common
#'   control in the secondary trial.
#' @param M Number of active treatment arms in the primary trial.
#' @param r Randomization ratio of each active arm to the common control
#'   in the primary trial.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is assumed.
#' @param L The interim adaptation look of the primary trial.
#' @param zL The z-test statistics at the interim adaptation look of
#'   the primary trial.
#' @param theta A vector of length \eqn{M} representing the assumed treatment
#'   effects for each active arm versus the common control. The global null
#'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
#'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
#' @param IMax Maximum information for any active arm versus the common
#'   control for the primary trial. Must be provided.
#' @param kMax The maximum number of stages of the primary trial.
#' @param informationRates The information rates of the primary trial.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param futilityStopping Indicators of whether futility stopping is
#'   allowed at each stage of the primary trial.
#' @param criticalValues The upper boundaries on the max z-test statistic
#'   scale for efficacy stopping for the primary trial. If missing, boundaries
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
#' @param futilityBounds The futility boundaries on the max-z statistic
#'   scale for the primary trial. Defaults to \code{rep(-8, kMax-1)}
#'   if left unspecified.
#' @param futilityCP The conditional power-based futility bounds for the
#'   primary trial.
#' @param futilityTheta The parameter value-based futility bounds for the
#'   primary trial.
#' @param spendingTime The error spending time of the primary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRates}.
#' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
#'   for trial adaptation.
#' @param MNew Number of active treatment arms in the secondary trial.
#' @param selected The indices of the selected active treatment arms for
#'   the secondary trial.
#' @param rNew Randomization ratio of each active arm to the common control
#'   in the secondary trial.
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
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param futilityBoundsInt The futility boundaries on the max-z statistic
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
#'   \code{"none"} for no early futility stopping.
#'   Defaults to \code{"none"}.
#' @param parameterBetaSpendingNew The parameter value of beta spending
#'   for the secondary trial. Corresponds to \eqn{\rho} for \code{"sfKD"},
#'   and \eqn{\gamma} for \code{"sfHSD"}.
#' @param spendingTimeNew The error spending time of the secondary trial.
#'   Defaults to missing, in which case it is assumed to be the same as
#'   \code{informationRatesNew}.
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
#' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
#' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
#'
#' @seealso \code{\link{adaptDesign_multiarm}}
#'
#' @examples
#'
#' getCP_multiarm(
#'   INew = 373 / 4, M = 2, r = 1, corr_known = FALSE,
#'   L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
#'   theta = c(-log(0.91), -log(0.78)),
#'   IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
#'   alpha = 0.025, typeAlphaSpending = "OF",
#'   MNew = 1, selected = 2, rNew = 1, nthreads = 1)
#'
#' @export
getCP_multiarm <- function(INew = NA_real_,
                           M = NA_integer_,
                           r = 1,
                           corr_known = TRUE,
                           L = NA_integer_,
                           zL = NA_real_,
                           theta = NA_real_,
                           IMax = NA_real_,
                           kMax = NA_integer_,
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
                           MNew = NA_integer_,
                           selected = NA_integer_,
                           rNew = 1,
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
                           nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getCP_multiarm_Rcpp(INew, M, r, corr_known, L, zL, theta, IMax,
                        kMax, informationRates, efficacyStopping,
                        futilityStopping, criticalValues,
                        alpha, typeAlphaSpending,
                        parameterAlphaSpending, userAlphaSpending,
                        futilityBounds, futilityCP, futilityTheta,
                        spendingTime, MullerSchafer, MNew, selected,
                        rNew, kNew, informationRatesNew,
                        efficacyStoppingNew, futilityStoppingNew,
                        typeAlphaSpendingNew,
                        parameterAlphaSpendingNew,
                        futilityBoundsInt, futilityCPInt,
                        futilityThetaInt, typeBetaSpendingNew,
                        parameterBetaSpendingNew, spendingTimeNew)
}


#' @title Confidence Interval After Trial Termination for Multi-Arm
#' Multi-Stage Design
#' @description Obtains the p-value, conservative point estimate, and
#' confidence interval after the end of a multi-arm multi-stage trial.
#'
#' @param M Number of active treatment arms.
#' @param r Randomization ratio of each active arm to the common control.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is used.
#' @param L The termination look.
#' @param zL The vector of z-test statistics at the termination look.
#' @param IMax Maximum information for any active arm versus the common
#'   control.
#' @param informationRates The information rates up to look \code{L}.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage up to look \code{L}.
#'   Defaults to \code{TRUE} if left unspecified.
#' @param criticalValues The matrix of by-level upper boundaries on the
#'   max z-test statistic scale for efficacy stopping up to look \code{L}.
#'   The first column is for level \code{M}, the second column is for
#'   level \code{M - 1}, and so on, with the last column for level 1.
#'   If left unspecified, the critical values will be computed based
#'   on the specified alpha spending function.
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
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' If \code{typeAlphaSpending} is \code{"OF"}, \code{"P"}, \code{"WT"}, or
#' \code{"none"}, then \code{informationRates}, \code{efficacyStopping},
#' and \code{spendingTime} must be of full length \code{kMax}, and
#' \code{informationRates} and \code{spendingTime} must end with 1.
#'
#' @return A data frame with the following components:
#'
#' * \code{level}: Number of elementary hypotheses considered for multiplicity.
#'
#' * \code{index}: The treatment arm with max Z among the active arms.
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
#' getCI_multiarm(
#'   L = 2, zL = c(2.075, 2.264),
#'   M = 2, r = 1, corr_known = FALSE,
#'   IMax = 300 / 4, informationRates = c(1/2, 1),
#'   alpha = 0.025, typeAlphaSpending = "sfOF",
#'   nthreads = 1)
#'
#' @export
getCI_multiarm <- function(M = NA_integer_,
                           r = 1,
                           corr_known = TRUE,
                           L = NA_integer_,
                           zL = NA_real_,
                           IMax = NA_real_,
                           informationRates = NA_real_,
                           efficacyStopping = NA_integer_,
                           criticalValues = NULL,
                           alpha = 0.025,
                           typeAlphaSpending = "sfOF",
                           parameterAlphaSpending = NA_real_,
                           spendingTime = NA_real_,
                           nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getCI_multiarm_Rcpp(M, r, corr_known, L, zL, IMax, informationRates,
                        efficacyStopping, criticalValues, alpha,
                        typeAlphaSpending, parameterAlphaSpending,
                        spendingTime)
}


#' @title Confidence Interval After Adaptation for Multi-Arm Multi-Stage
#' Design
#' @description Obtains the p-value, conservative point estimate, and
#' confidence interval after the end of an adaptive multi-arm multi-stage trial.
#'
#' @param M Number of active treatment arms in the primary trial.
#' @param r Randomization ratio of each active arm to the common control
#'   in the primary trial.
#' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
#'   statistics is derived from the randomization ratio \eqn{r}
#'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
#'   0 is assumed.
#' @param L The interim adaptation look of the primary trial.
#' @param zL The z-test statistics at the interim adaptation look of
#'   the primary trial.
#' @param IMax Maximum information for any active arm versus the common
#'   control for the primary trial. Must be provided.
#' @param kMax The maximum number of stages of the primary trial.
#' @param informationRates The information rates of the primary trial.
#' @param efficacyStopping Indicators of whether efficacy stopping is
#'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
#'   if left unspecified.
#' @param criticalValues The matrix of by-level upper boundaries on the
#'   max z-test statistic scale for efficacy stopping up to look \code{L}
#'   for the primary trial.
#'   The first column is for level \code{M}, the second column is for
#'   level \code{M - 1}, and so on, with the last column for level 1.
#'   If left unspecified, the critical values will be computed based
#'   on the specified alpha spending function.
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
#' @param MNew The number of active treatment arms in the secondary trial.
#' @param selected The indices of the selected treatment arms for the
#'   secondary trial among the \code{M} active arms in the primary trial.
#' @param rNew The randomization ratio of each active arm to the common control
#'  in the secondary trial.
#' @param Lc The termination look of the integrated trial.
#' @param zLc The z-test statistics at the termination look of the
#'   integrated trial.
#' @param INew The maximum information for any active arm versus the common
#'   control in the secondary trial.
#' @param informationRatesNew The spacing of looks of the secondary trial.
#' @param efficacyStoppingNew The indicators of whether efficacy stopping is
#'   allowed at each look of the secondary trial.
#'   Defaults to \code{TRUE} if left unspecified.
#' @param typeAlphaSpendingNew The type of alpha spending for the secondary
#'   trial. One of the following:
#'   \code{"OF"} for O'Brien-Fleming boundaries,
#'   \code{"sfOF"} for O'Brien-Fleming type spending function,
#'   \code{"sfP"} for Pocock type spending function,
#'   \code{"sfKD"} for Kim & DeMets spending function,
#'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
#'   \code{"none"} for no early efficacy stopping.
#'   Defaults to \code{"sfOF"}.
#' @param parameterAlphaSpendingNew The parameter value of alpha spending
#'   for the secondary trial. Corresponds to
#'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
#' @param spendingTimeNew The error spending time of the secondary trial.
#'   Defaults to missing, in which case, it is
#'   the same as \code{informationRatesNew}.
#' @param nthreads The number of threads to use (0 leaves the
#'   RcppParallel setting unchanged).
#'
#' @details
#' If typeAlphaSpendingNew is \code{"OF"} or \code{"none"}, then
#' \code{informationRatesNew}, \code{efficacyStoppingNew}, and
#' \code{spendingTimeNew} must be of full length \code{kNew}, and
#' \code{informationRatesNew} and \code{spendingTimeNew} must end with 1.
#'
#' @return A data frame with the following variables:
#'
#' * \code{level}: Number of elementary hypotheses considered for multiplicity.
#'
#' * \code{index}: The treatment arm with max Z among the active arms.
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
#' getADCI_multiarm(
#'   M = 2, r = 1, corr_known = FALSE, L = 1, zL = c(2.075, 2.264),
#'   IMax = 300 / 4, kMax = 2, informationRates = c(0.5, 1),
#'   alpha = 0.025, typeAlphaSpending = "sfOF",
#'   MNew = 1, selected = 2, rNew = 1,
#'   Lc = 2, zLc = 1.667, INew = 374 / 4, nthreads = 1)
#'
#' @export
getADCI_multiarm <- function(M = NA_integer_,
                             r = 1,
                             corr_known = TRUE,
                             L = NA_integer_,
                             zL = NA_real_,
                             IMax = NA_real_,
                             kMax = NA_integer_,
                             informationRates = NA_real_,
                             efficacyStopping = NA_integer_,
                             criticalValues = NULL,
                             alpha = 0.25,
                             typeAlphaSpending = "sfOF",
                             parameterAlphaSpending = NA_real_,
                             spendingTime = NA_real_,
                             MullerSchafer = FALSE,
                             MNew = NA_integer_,
                             selected = NA_integer_,
                             rNew = 1,
                             Lc = NA_integer_,
                             zLc = NA_real_,
                             INew = NA_real_,
                             informationRatesNew = NA_real_,
                             efficacyStoppingNew = NA_integer_,
                             typeAlphaSpendingNew = "sfOF",
                             parameterAlphaSpendingNew = NA_real_,
                             spendingTimeNew = NA_real_,
                             nthreads = 0) {
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }
    getADCI_multiarm_Rcpp(M, r, corr_known, L, zL, IMax, kMax,
                          informationRates, efficacyStopping,
                          criticalValues, alpha, typeAlphaSpending,
                          parameterAlphaSpending, spendingTime,
                          MullerSchafer, MNew, selected, rNew,
                          Lc, zLc, INew, informationRatesNew,
                          efficacyStoppingNew, typeAlphaSpendingNew,
                          parameterAlphaSpendingNew,
                          spendingTimeNew)
}
