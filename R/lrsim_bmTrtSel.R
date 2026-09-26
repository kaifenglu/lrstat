#' @title Simulation of a seamless phase II/III design with treatment
#'   selection based on a short-term endpoint and toxicities
#' @description Simulates a two-stage seamless phase II/III trial in which
#'   several doses are compared with a common control. At the end of phase II
#'   a single dose is carried forward based on the posterior benefit-risk
#'   tradeoff of a binary short-term biomarker endpoint and a binary toxicity
#'   endpoint. The confirmatory phase III analysis is performed on a
#'   time-to-event long-term endpoint linked to the biomarker through the
#'   copula correlation, and the type I error rate is protected by a closed
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
#' @param corrEfficacyToxicity The correlation between the latent normal
#'   variable for the toxicity endpoint and \eqn{z_B}, where \eqn{z_B} is the
#'   latent normal variable that determines the short-term biomarker endpoint.
#'   Use 0 for independent biomarker and toxicity endpoints.
#' @param corrEfficacyTTE The correlation between the latent normal variable for
#'   the long-term time-to-event endpoint and \eqn{z_B}. Use 0 for an
#'   independent biomarker and long-term endpoint.
#' @param hazardRateControl The overall hazard rate for the long-term
#'   endpoint in the control arm. It must be positive.
#' @param hazardRatioTreatments A vector of length \code{M} giving the hazard
#'   ratio for each treatment relative to the control arm. The hazard rate in
#'   treatment arm \code{k} is
#'   \code{hazardRateControl * hazardRatioTreatments[k]}; values below 1
#'   indicate effective treatments.
#' @param totalNumberOfEvents The target total number of long-term endpoint
#'   events in the selected dose and control arms, including subjects enrolled
#'   in phases II and III. It must be a positive integer no greater than the
#'   combined sample size of those two arms at
#'   \code{phase3SampleSizePerArmMin}. When provided, the final analysis is
#'   event-driven and \code{studyDurationPhase3} is ignored.
#' @param studyDurationPhase3 The duration of phase III, measured from the
#'   opening of phase III enrollment to the final analysis. It must be provided
#'   when \code{totalNumberOfEvents} is missing.
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
#' @param useUniformPrior Whether to use the uniform Beta(1,1) prior (the
#'   default) or the Jeffreys Beta(0.5,0.5) prior for the beta-binomial
#'   posterior used in dose selection.
#' @param methods A character vector naming the testing procedures to evaluate
#'   for the confirmatory analysis. Any subset of \code{"ctbonferroni"},
#'   \code{"ctdunnett"}, \code{"ctsimes"}, \code{"ctpooled"}, \code{"cer"},
#'   \code{"tsssd.k"},
#'   \code{"tsssd.uk"}, \code{"tsssd.k.rank"}, \code{"tsssd.uk.rank"},
#'   \code{"tsssd.k.ce"}, \code{"tsssd.uk.ce"},
#'   \code{"tsssd.k.rank.ce"}, \code{"tsssd.uk.rank.ce"},
#'   \code{"bm.rank"}, \code{"pe.rank"},
#'   \code{"naive"}, and \code{"ph3only"}. Restricting the set skips the
#'   corresponding computation entirely, which matters because the methods
#'   differ by orders of magnitude in cost.
#' @param accrualRatePhase2 The accrual rate per arm during phase II. Arrival
#'   times follow a homogeneous Poisson process.
#' @param accrualRatePhase3 The accrual rate per arm during phase III.
#' @param followupTimePhase2 The follow-up time after the last phase II
#'   enrollment across all arms. Phase III enrollment opens at that point. Use
#'   0 when dose selection occurs immediately after the last phase II
#'   enrollment.
#' @param maxNumberOfIterations The number of simulated trials.
#' @param maxNumberOfRawDatasets The number of initial simulation iterations
#'   for which to retain subject-level raw data. Set to 0, the default, to skip
#'   raw-data retention.
#' @param seed The seed for the random number generator.
#' @param nthreads The number of threads to use. The default, 0, leaves the
#'   \code{RcppParallel} setting unchanged.
#'
#' @return A list of operating characteristics. Let \code{ngrid} denote
#'   \code{phase3SampleSizePerArmMax - phase3SampleSizePerArmMin + 1}, the
#'   number of stage 2 sample sizes examined, and let \code{M} denote the
#'   number of doses. The list contains
#'
#' * \code{n1}, \code{n2}, \code{totalNumberOfEvents},
#'   \code{studyDurationPhase3}, \code{numberOfIterations}, \code{trueOBD}:
#'   The design inputs echoed back, with \code{n2} the vector of stage 2 sample
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
#'   - \code{disjunctive.power}: A vector of length \code{ngrid} giving the
#'     probability of rejecting at least one true non-null hypothesis. It is
#'     missing when all treatment hazard ratios equal 1.
#'
#'   - \code{conjunctive.power}: A vector of length \code{ngrid} giving the
#'     probability of rejecting all true non-null hypotheses. It is missing
#'     when all treatment hazard ratios equal 1.
#'
#'   - \code{fwer}: A vector of length \code{ngrid} giving the probability
#'     of rejecting any true null hypothesis.
#'
#'   The method names are
#'   \code{ctbonferroni}, \code{ctdunnett}, \code{ctsimes}, and
#'   \code{ctpooled} for the closed testing procedure with the inverse normal
#'   combination of stage 1 and stage 2 p-values, using the Bonferroni,
#'   Dunnett, Simes, and pooled log-rank local tests respectively; \code{cer}
#'   for the conditional error rate method;
#'   \code{tsssd.k} and \code{tsssd.uk} for the original two-stage seamless
#'   design boundaries with known and unknown correlation; \code{tsssd.k.rank}
#'   and \code{tsssd.uk.rank} for rank-based boundaries based on the effective
#'   number of less efficacious doses; and \code{tsssd.k.ce},
#'   \code{tsssd.uk.ce}, \code{tsssd.k.rank.ce}, and \code{tsssd.uk.rank.ce}
#'   for conditional-error updates of the original and rank-based boundaries.
#'   These updates start from nominal boundaries based
#'   on \code{n1/(n1+n2)} and then use the observed stage 1 z-statistic and
#'   observed information fraction;
#'   \code{bm.rank} for the Biomarker rank-based Dunnett adjustment for unknown
#'   biomarker-efficacy correlation in Wang et al. formula (3) for stage 1,
#'   along with the p-value combination test at the end of stage 2;
#'   \code{pe.rank} for the primary endpoint rank-based Dunnett adjustment for
#'   stage 1, along with the p-value combination test at the end of stage 2;
#'   \code{naive} for the
#'   unadjusted log-rank test on the combined stage 1 and stage 2 data; and
#'   \code{ph3only} for the unadjusted log-rank test on the stage 2 data only.
#'   Because dose selection uses stage 1 data only, \code{ph3only} is based on
#'   data independent of the selection and still controls the familywise error
#'   rate, at the cost of discarding the stage 1 information. \code{naive}
#'   reuses the selection data and is anticonservative; it is reported for
#'   reference.
#'
#' * \code{sumdataBIN}: One row per iteration and phase II arm, containing
#'   \code{iterationNumber}, \code{treatmentGroup}, \code{responses},
#'   \code{toxicities}, \code{zBiomarker}, and \code{selected}. Treatment
#'   group 0 is control; its toxicity count and biomarker Z statistic are
#'   missing. These rows reproduce the dose-selection operating
#'   characteristics.
#'
#' * \code{sumdataTTE}: One row per iteration and phase III sample size,
#'   containing \code{iterationNumber}, \code{selectedDose},
#'   \code{phase3SampleSize}, stage 1, stage 2, and total event counts,
#'   stage 1, stage 2, and cumulative log-rank Z statistics for the selected
#'   dose, as well as stage 1 log-rank Z statistics for all individual and
#'   pooled dose groups (e.g., \code{stage1LogRankZ1}, \code{stage1LogRankZ2},
#'   \code{stage1LogRankZ12}). It also has one \code{reject.<method>}
#'   indicator column for every requested method. These rows reproduce
#'   \code{ave.event} and the operating characteristics in \code{byMethod}.
#'
#' * \code{rawdataBIN} (present when \code{maxNumberOfRawDatasets} is
#'   positive): Subject-level phase II binary data with \code{iterationNumber},
#'   \code{subjectId}, \code{treatmentGroup} (0 for control and 1 through
#'   \code{M} for doses), \code{response}, and \code{toxicity}. Toxicity is
#'   missing for the control arm because it is not simulated or used for dose
#'   selection.
#'
#' * \code{rawdataTTE} (present when \code{maxNumberOfRawDatasets} is
#'   positive): Subject-level time-to-event data with \code{iterationNumber},
#'   \code{subjectId}, \code{phase}, \code{treatmentGroup},
#'   \code{response}, \code{arrivalTime}, \code{survivalTime},
#'   \code{timeUnderObservation}, and \code{event}. Phase 2 includes all
#'   arms; phase 3 includes only the control arm and the selected dose.
#'
#' @details
#' Data generation uses a copula-based approach. For each subject, a base
#' random variable \eqn{z_B \sim N(0,1)} drives the short-term biomarker
#' endpoint, \code{shortv}. Toxicity is generated from
#' \eqn{z_S = \rho_{tox}z_B + \sqrt{1-\rho_{tox}^2}z_{S,indep}} and the
#' toxicity indicator \code{toxv} is obtained by thresholding \eqn{z_S} at
#' the specified toxicity probability. The long-term endpoint is generated
#' from \eqn{z_E = \rho_{eff}z_B + \sqrt{1-\rho_{eff}^2}z_{E,indep}} as
#' \eqn{v = -\log(\Phi(z_E))/\lambda}, where \eqn{\lambda =
#' hazardRateControl} in the control arm and
#' \eqn{\lambda = hazardRateControl \times hazardRatioTreatments[k]}
#' in treatment arm \eqn{k}. Dose selection uses a beta-binomial model with
#' independent priors, either uniform Beta(1,1) or Jeffreys Beta(0.5,0.5)
#' depending on \code{useUniformPrior}. A dose enters the acceptable set when
#' the posterior probability that its response rate exceeds that of the
#' control is above \code{efficacyThreshold} and the posterior probability
#' that its toxicity rate is below \code{toxicityUpperLimit} is above
#' \code{safetyThreshold}. Among the acceptable doses, the one maximizing the
#' posterior mean benefit-risk tradeoff is selected. When both thresholds are
#' 0, all doses are acceptable.
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
#' Ping Gao and Yingqiu Li. Adaptive two-stage seamless sequential design for
#' clinical trials. Journal of Biopharmaceutical Statistics. 2025, 35(4),
#' 565-587.
#'
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#'
#' # Overall control-arm hazard rate for the long-term endpoint
#' hazardRateControl <- log(2) / 15.5
#'
#' # response rates of the two doses under investigation
#' pe <- c(0.6, 0.5)
#'
#' # Hazard ratio versus control for each treatment
#' hazardRatioTreatments <- c(0.65, 0.70)
#'
#' sim <- lrsim_bmTrtSel(
#'   phase2SampleSizePerArm = 50,
#'   phase3SampleSizePerArmMin = 113,
#'   phase3SampleSizePerArmMax = 118,
#'   responseProbControl = 0.4,
#'   responseProbTreatments = pe,
#'   toxicityProbTreatments = c(0, 0),
#'   corrEfficacyToxicity = 0,
#'   corrEfficacyTTE = 0.43,
#'   hazardRateControl = hazardRateControl,
#'   hazardRatioTreatments = hazardRatioTreatments,
#'   studyDurationPhase3 = 42.1,
#'   toxicityWeight = 0,
#'   toxicityUpperLimit = 1,
#'   efficacyThreshold = 0,
#'   safetyThreshold = 0,
#'   methods = c("ctbonferroni", "ctdunnett", "ctsimes", "ctpooled",
#'               "cer", "naive", "ph3only"),
#'   accrualRatePhase2 = 3,
#'   accrualRatePhase3 = 6,
#'   followupTimePhase2 = 6,
#'   maxNumberOfIterations = 100,
#'   seed = 314159,
#'   nthreads = 1)
#'
#' sim$pcs
#' sim$byMethod$ctdunnett$gpower
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
    corrEfficacyTTE = 0,
    hazardRateControl = NA_real_,
    hazardRatioTreatments = NA_real_,
    totalNumberOfEvents = NA_integer_,
    studyDurationPhase3 = NA_real_,
    toxicityWeight = NA_real_,
    toxicityUpperLimit = NA_real_,
    efficacyThreshold = 0,
    safetyThreshold = 0,
    useUniformPrior = TRUE,
    methods = c("ctbonferroni", "ctdunnett", "ctsimes", "ctpooled", "cer",
                "tsssd.k", "tsssd.uk", "tsssd.k.rank", "tsssd.uk.rank",
                "tsssd.k.ce", "tsssd.uk.ce",
                "tsssd.k.rank.ce", "tsssd.uk.rank.ce",
                "bm.rank", "pe.rank",
                "naive", "ph3only"),
    accrualRatePhase2 = NA_real_,
    accrualRatePhase3 = NA_real_,
    followupTimePhase2 = 0,
    maxNumberOfIterations = 1000,
    maxNumberOfRawDatasets = 0,
    seed = 0,
    nthreads = 0) {

  # Respect user-requested number of threads (best effort)
  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  methods <- tolower(methods)

  lrsim_bmTrtSel_Rcpp(
    phase2SampleSizePerArm = phase2SampleSizePerArm,
    phase3SampleSizePerArmMin = phase3SampleSizePerArmMin,
    phase3SampleSizePerArmMax = phase3SampleSizePerArmMax,
    responseProbControl = responseProbControl,
    responseProbTreatments = responseProbTreatments,
    toxicityProbTreatments = toxicityProbTreatments,
    corrEfficacyToxicity = corrEfficacyToxicity,
    corrEfficacyTTE = corrEfficacyTTE,
    hazardRateControl = hazardRateControl,
    hazardRatioTreatments = hazardRatioTreatments,
    totalNumberOfEvents = totalNumberOfEvents,
    studyDurationPhase3 = studyDurationPhase3,
    toxicityWeight = toxicityWeight,
    toxicityUpperLimit = toxicityUpperLimit,
    efficacyThreshold = efficacyThreshold,
    safetyThreshold = safetyThreshold,
    useUniformPrior = useUniformPrior,
    methods = methods,
    accrualRatePhase2 = accrualRatePhase2,
    accrualRatePhase3 = accrualRatePhase3,
    followupTimePhase2 = followupTimePhase2,
    maxNumberOfIterations = maxNumberOfIterations,
    maxNumberOfRawDatasets = maxNumberOfRawDatasets,
    seed = seed)
}
