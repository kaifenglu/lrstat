#' Simulate the Carreras adaptive seamless design with PK-guided selection
#'
#' Simulates the oncology case study of Carreras, Gutjahr, and Brannath
#' (2015). Two experimental regimens and a control are randomized in stage 1,
#' one regimen is selected at an interim analysis, and the selected regimen
#' and control continue in stage 2. Overall survival is tested using the
#' follow-up-wise, patient-wise, conservative Dunnett, stage-2-only, and
#' conventional three-arm procedures considered in the article.
#'
#' @param medianSurvival Median overall survival in months for the high-dose,
#'   low-dose, and control arms, in that order.
#' @param meanExposure Mean AUC for the high-dose and low-dose arms.
#' @param exposureSD Common standard deviation of AUC on its original scale.
#' @param corrSurvivalExposure Gaussian copula correlation between survival
#'   time and log exposure.
#' @param interimPatients Number of patients in the stage-1 cohort. The interim
#'   analysis occurs after the last of these patients has at least
#'   `interimFollowup` months of follow-up. Accrual continues in all three arms
#'   until that analysis; these additional recruits belong to stage 2 and are
#'   excluded from interim OS statistics, but their available exposure data
#'   contribute to exposure-based treatment selection.
#' @param interimFollowup Minimum follow-up in months at the interim analysis.
#' @param totalPatients Total sample size across both stages.
#' @param finalAnalysisTime Calendar time in months from first enrollment to
#'   the final analysis.
#' @param selectionRule Regimen-selection rule: `"exposure"` selects high dose
#'   when its observed mean AUC is at least 1.5 times the low-dose mean;
#'   `"os"` selects the regimen with the better interim log-rank statistic;
#'   `"perfect"` selects the regimen with the smaller final hazard-ratio
#'   estimate across the two counterfactual adaptive continuations.
#' @param maxNumberOfIterations Number of simulated trials.
#' @param seed Random-number seed.
#' @param nthreads Number of simulation threads. Zero leaves the current
#'   RcppParallel setting unchanged.
#'
#' @return An object of class `lrsim_pkTrtSel`. It contains selection
#'   probabilities, selection bias in the hazard-ratio estimate, estimated
#'   combination-test weights, rejection probabilities by method, and a
#'   trial-level `sumdata` data frame.
#'
#' @details
#' Accrual follows the article's piecewise rates of 7 patients/month for the
#' first 4 months, 15 patients/month for the next 4 months, and 22
#' patients/month thereafter. Stage 1 uses 1:1:1 allocation. Stage-2 patients
#' accrued before treatment selection also use 1:1:1 allocation; subsequent
#' patients use 1:1 allocation to control and the selected regimen. Exposure is lognormal,
#' parameterized by the supplied arithmetic mean and standard deviation, and
#' is joined to exponential survival through a Gaussian copula. Combination
#' weights are estimated from the average simulated event counts as described
#' in the article.
#'
#' @references
#' Carreras M, Gutjahr G, Brannath W. Adaptive seamless designs with interim
#' treatment selection: a case study in oncology. *Statistics in Medicine*.
#' 2015;34:1317-1333. \doi{10.1002/sim.6407}.
#'
#' @examples
#' sim <- lrsim_pkTrtSel(maxNumberOfIterations = 20, seed = 314159,
#'                       nthreads = 1)
#' sim$probabilitySelectHigh
#' sim$byMethod$dunnett$rejectionProbability
#'
#' @export
lrsim_pkTrtSel <- function(
    medianSurvival = c(9, 7.5, 6),
    meanExposure = c(300, 200),
    exposureSD = 120,
    corrSurvivalExposure = 0.5,
    interimPatients = 100L,
    interimFollowup = 3,
    totalPatients = 410L,
    finalAnalysisTime = 29,
    selectionRule = c("exposure", "os", "perfect"),
    maxNumberOfIterations = 1000L,
    seed = 0L,
    nthreads = 0L) {
  selectionRule <- match.arg(tolower(selectionRule),
                             c("exposure", "os", "perfect"))
  if (nthreads > 0L) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  lrsim_pkTrtSel_Rcpp(
    medianSurvival = medianSurvival,
    meanExposure = meanExposure,
    exposureSD = exposureSD,
    corrSurvivalExposure = corrSurvivalExposure,
    interimPatients = interimPatients,
    interimFollowup = interimFollowup,
    totalPatients = totalPatients,
    finalAnalysisTime = finalAnalysisTime,
    selectionRule = selectionRule,
    maxNumberOfIterations = maxNumberOfIterations,
    seed = seed
  )
}