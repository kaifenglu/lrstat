#' @title Simulation of a seamless phase II/III design with treatment
#'   selection based on a short-term endpoint
#' @description Simulates a two-stage seamless phase II/III trial in which
#'   several doses are compared with a common control. At the end of phase II
#'   a single dose is carried forward based on the posterior benefit-risk
#'   tradeoff of a binary short-term efficacy endpoint and a binary toxicity
#'   endpoint. The confirmatory phase III analysis is performed on a
#'   time-to-event long-term endpoint whose hazard depends on the short-term
#'   response status, and the type I error rate is protected by a closed
#'   testing procedure combined across the two stages.
#'
#' @param phase2SampleSizePerArm The number of subjects per arm enrolled in
#'   phase II (stage 1).
#' @param phase3SampleSizePerArmMin The smallest number of subjects per arm
#'   enrolled in phase III (stage 2). Operating characteristics are reported
#'   for every stage 2 sample size from this value to
#'   \code{phase3SampleSizePerArmMax}.
#' @param phase3SampleSizePerArmMax The largest number of subjects per arm
#'   enrolled in phase III (stage 2).
#' @param responseProbControl The probability of a short-term response in the
#'   control arm.
#' @param responseProbTreatments A vector of length \code{M} giving the
#'   probability of a short-term response for each of the \code{M} doses under
#'   investigation. Its length determines \code{M}.
#' @param toxicityProbTreatments A vector of length \code{M} giving the
#'   probability of toxicity for each dose under investigation.
#' @param corrEfficacyToxicity The correlation between the bivariate latent
#'   normal variables used to generate the binary efficacy and toxicity
#'   endpoints. This is the correlation on the latent scale, not the
#'   correlation of the observed binary endpoints. Use 0 for independent
#'   endpoints.
#' @param hazardRateControl A vector of length 2 giving the hazard rate of the
#'   long-term endpoint in the control arm for short-term nonresponders and
#'   responders, respectively.
#' @param hazardRateTreatments An \code{M} by 2 matrix giving the hazard rate
#'   of the long-term endpoint for each dose, with the first column for
#'   short-term nonresponders and the second column for responders.
#' @param studyDurationPhase3 The duration of phase III, measured from the
#'   first phase III enrollment to the final analysis.
#' @param toxicityWeight The weight placed on the posterior mean toxicity rate
#'   in the benefit-risk tradeoff used for dose selection. Use 0 to select on
#'   efficacy alone.
#' @param toxicityUpperLimit The prespecified upper limit for the toxicity
#'   rate used in the safety criterion. Use 1 when the safety criterion is not
#'   applied.
#' @param efficacyThreshold The threshold for the posterior probability that a
#'   dose is superior to the control in short-term response. Use 0 when the
#'   efficacy criterion is not applied.
#' @param safetyThreshold The threshold for the posterior probability that the
#'   toxicity rate of a dose is below \code{toxicityUpperLimit}. Use 0 when
#'   the safety criterion is not applied.
#' @param methods A character vector naming the testing procedures to evaluate
#'   for the confirmatory analysis. Any subset of \code{"ctdunnett"},
#'   \code{"ctsimes"}, \code{"ctpooled"}, \code{"cer"}, \code{"TSSSP.k"},
#'   \code{"TSSSP.uk"}, \code{"naive"}, and \code{"ph3only"}. Restricting the
#'   set skips the corresponding computation entirely, which matters because
#'   the methods differ by orders of magnitude in cost.
#' @param accrualRatePhase2 The accrual rate per arm during phase II. Arrival
#'   times follow a homogeneous Poisson process.
#' @param accrualRatePhase3 The accrual rate per arm during phase III.
#' @param followupTimePhase2 The follow-up time after the last phase II
#'   enrollment across all arms. Phase III enrollment opens at that point. Use
#'   0 when dose selection occurs immediately after the last phase II
#'   enrollment.
#' @param maxNumberOfIterations The number of simulated trials.
#' @param seed The seed for the random number generator.
#' @param nthreads The number of threads to use. The default, 0, leaves the
#'   \code{RcppParallel} setting unchanged.
#'
#' @return A list of operating characteristics. Let \code{ngrid} denote
#'   \code{phase3SampleSizePerArmMax - phase3SampleSizePerArmMin + 1}, the
#'   number of stage 2 sample sizes examined, and let \code{M} denote the
#'   number of doses. The list contains
#'
#' * \code{n1}, \code{n2}, \code{numberOfIterations}, \code{trueOBD}: The
#'   design inputs echoed back, with \code{n2} the vector of stage 2 sample
#'   sizes examined.
#'
#' * \code{selectionProb}: A vector of length \code{M} giving the probability
#'   that each dose is selected at the end of phase II.
#'
#' * \code{pcs}: The percentage of simulated trials selecting the dose with
#'   the largest true benefit-risk tradeoff \code{responseProbTreatments -
#'   toxicityWeight * toxicityProbTreatments}.
#'
#' * \code{ave.event}: An \code{ngrid} by 3 matrix of the average number of
#'   events in the selected dose and the control arm combined, in stage 1,
#'   stage 2, and overall.
#'
#' * \code{methods}: The testing procedures evaluated, in canonical order.
#'
#' * \code{byMethod}: A named list with one element per evaluated method,
#'   each containing
#'
#'   - \code{gpower}: A vector of length \code{ngrid} giving the generalized
#'     power, the probability of both selecting the true best dose and
#'     rejecting its null hypothesis.
#'
#'   - \code{prob.rej.each}: An \code{ngrid} by \code{M} matrix of the
#'     probability of rejecting the null hypothesis for each dose conditional
#'     on that dose being selected.
#'
#'   - \code{prob.rej.any}: A vector of length \code{ngrid} giving the
#'     probability of rejecting any null hypothesis.
#'
#'   - \code{elapsedTime}: The total thread time in seconds attributed to the
#'     method, summed over all simulated trials and stage 2 sample sizes. Work
#'     shared across methods, in particular the log-rank tests and the data
#'     generation, is excluded, so these are incremental costs rather than a
#'     decomposition of the total run time.
#'
#'   The method names are
#'   \code{ctdunnett}, \code{ctsimes}, and \code{ctpooled} for the closed
#'   testing procedure with the inverse normal combination of stage 1 and
#'   stage 2 p-values, using the Dunnett, Simes, and pooled log-rank local
#'   tests respectively; \code{cer} for the conditional error rate method;
#'   \code{TSSSP.k} and \code{TSSSP.uk} for the two-stage seamless design
#'   boundaries with known and unknown correlation; \code{naive} for the
#'   unadjusted log-rank test on the combined stage 1 and stage 2 data; and
#'   \code{ph3only} for the unadjusted log-rank test on the stage 2 data only.
#'   Because dose selection uses stage 1 data only, \code{ph3only} is based on
#'   data independent of the selection and still controls the familywise error
#'   rate, at the cost of discarding the stage 1 information. \code{naive}
#'   reuses the selection data and is anticonservative; it is reported for
#'   reference.
#'
#' @details
#' For each subject on a treatment arm, a bivariate latent normal vector
#' \eqn{(z_T, z_E)} with mean zero, unit variances, and correlation
#' \code{corrEfficacyToxicity} is drawn. The binary toxicity and efficacy
#' endpoints are obtained as \eqn{Y_T = I\{z_T \le \Phi^{-1}(p_T(d))\}} and
#' \eqn{Y_E = I\{z_E \le \Phi^{-1}(p_E(d))\}}, so that the marginal
#' probabilities are \eqn{p_T(d)} and \eqn{p_E(d)} while the two endpoints are
#' correlated. The long-term endpoint is exponential with a rate determined by
#' the realized short-term response status.
#'
#' Dose selection uses a beta-binomial model with independent uniform priors.
#' A dose enters the acceptable set when the posterior probability that its
#' response rate exceeds that of the control is above \code{efficacyThreshold}
#' and the posterior probability that its toxicity rate is below
#' \code{toxicityUpperLimit} is above \code{safetyThreshold}. Among the
#' acceptable doses, the one maximizing the posterior mean benefit-risk
#' tradeoff is selected. When both thresholds are 0, all doses are acceptable.
#'
#' Phase III enrollment opens \code{followupTimePhase2} after the last phase
#' II enrollment across all arms, and the final analysis occurs
#' \code{studyDurationPhase3} later. Subjects whose long-term endpoint has not
#' occurred by then are censored at the analysis time.
#'
#' @references
#' Liyun Jiang and Ying Yuan. Seamless phase II/III design: a useful strategy
#' to reduce the sample size for dose optimization. Journal of the National
#' Cancer Institute. 2023, 115(9):1092-1098.
#'
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # hazard rates in the nonresponse and response groups of the control arm
#' the0 <- c(log(2)/12, log(2)/24)
#'
#' # response rates of the two doses under investigation
#' pe <- c(0.6, 0.5)
#'
#' # hazard ratio versus control within each response group
#' hr <- rbind(c(0.75, 0.75), c(0.75, 0.75))
#' the1 <- t(sapply(1:2, function(k) hr[k,]*the0))
#'
#' sim <- lrsim_bmTrtSel(
#'   phase2SampleSizePerArm = 50,
#'   phase3SampleSizePerArmMin = 113,
#'   phase3SampleSizePerArmMax = 118,
#'   responseProbControl = 0.4,
#'   responseProbTreatments = pe,
#'   toxicityProbTreatments = c(0, 0),
#'   corrEfficacyToxicity = 0,
#'   hazardRateControl = the0,
#'   hazardRateTreatments = the1,
#'   studyDurationPhase3 = 42.1,
#'   toxicityWeight = 0,
#'   toxicityUpperLimit = 1,
#'   efficacyThreshold = 0,
#'   safetyThreshold = 0,
#'   accrualRatePhase2 = 3,
#'   accrualRatePhase3 = 6,
#'   followupTimePhase2 = 6,
#'   maxNumberOfIterations = 100,
#'   seed = 314159)
#'
#' sim$pcs
#' sim$byMethod$ctdunnett$gpower
#'
#' # timing comparison across methods
#' sapply(sim$byMethod, function(m) m$elapsedTime)
#'
#' @export
lrsim_bmTrtSel <- function(
    phase2SampleSizePerArm = NA_integer_,
    phase3SampleSizePerArmMin = NA_integer_,
    phase3SampleSizePerArmMax = NA_integer_,
    responseProbControl = NA_real_,
    responseProbTreatments = NA_real_,
    toxicityProbTreatments = NA_real_,
    corrEfficacyToxicity = 0,
    hazardRateControl = NA_real_,
    hazardRateTreatments = matrix(),
    studyDurationPhase3 = NA_real_,
    toxicityWeight = NA_real_,
    toxicityUpperLimit = NA_real_,
    efficacyThreshold = 0,
    safetyThreshold = 0,
    methods = c("ctdunnett", "ctsimes", "ctpooled", "cer",
                "TSSSP.k", "TSSSP.uk", "naive", "ph3only"),
    accrualRatePhase2 = NA_real_,
    accrualRatePhase3 = NA_real_,
    followupTimePhase2 = 0,
    maxNumberOfIterations = 1000,
    seed = 0,
    nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsim_bmTrtSel_Rcpp(
    phase2SampleSizePerArm = phase2SampleSizePerArm,
    phase3SampleSizePerArmMin = phase3SampleSizePerArmMin,
    phase3SampleSizePerArmMax = phase3SampleSizePerArmMax,
    responseProbControl = responseProbControl,
    responseProbTreatments = responseProbTreatments,
    toxicityProbTreatments = toxicityProbTreatments,
    corrEfficacyToxicity = corrEfficacyToxicity,
    hazardRateControl = hazardRateControl,
    hazardRateTreatments = hazardRateTreatments,
    studyDurationPhase3 = studyDurationPhase3,
    toxicityWeight = toxicityWeight,
    toxicityUpperLimit = toxicityUpperLimit,
    efficacyThreshold = efficacyThreshold,
    safetyThreshold = safetyThreshold,
    methods = methods,
    accrualRatePhase2 = accrualRatePhase2,
    accrualRatePhase3 = accrualRatePhase3,
    followupTimePhase2 = followupTimePhase2,
    maxNumberOfIterations = maxNumberOfIterations,
    seed = seed)
}
