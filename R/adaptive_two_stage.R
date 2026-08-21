
#' @title Compute stage-wise local p-values
#' @description Compute local p-values for the intersection hypotheses in a
#'   two-stage adaptive multiple testing procedure using a p-value combination
#'   method.
#' @param stg2_p Stage 2 p-values for the elementary hypotheses.
#' @param wgtmat Weight matrix for the stage 1 intersection hypotheses. If
#'   \code{NULL}, equal weights are assigned within each intersection
#'   hypothesis.
#' @param family Family matrix indicating which hypotheses belong to which
#'   families. The correlation is known only for hypotheses belonging to the
#'   same family. Defaults to one family containing all elementary hypotheses.
#' @param corr Correlation matrix for the test statistics. If \code{NULL},
#'   within-family correlations are 0.5 and between-family correlations are
#'   missing.
#' @param stg1_inthyp_nr Indices of the stage 1 intersection hypotheses that
#'   were not rejected.
#' @param stg2_elemhyp Indices of the elementary hypotheses tested at stage 2.
#' @param stg2_wgtmat Weight matrix for the stage 2 intersection hypotheses.
#'   If \code{NULL}, equal weights are assigned within each intersection
#'   hypothesis.
#' @param test P-value combination method. It can start with \code{"bon"},
#'   \code{"sim"}, or \code{"dun"}; the default is \code{"dunnett"}.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list containing \code{inthyp_idx}, the indices of the stage 1
#'   intersection hypotheses, \code{inthyp}, their indicator matrix, and
#'   \code{pinter}, their local p-values.
#'
#' @details
#' Despite the stage-oriented parameter names, this function can also be used
#' to generate stage 1 local p-values. In that case, provide the complete set
#' of intersection hypotheses in \code{stg1_inthyp_nr}, all elementary
#' hypotheses in \code{stg2_elemhyp}, and the corresponding stage 1 p-values
#' and weight matrix. The example illustrates this use.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' fadjp(c(0.00045, 0.0952, 0.0225, 0.1104),
#'       wgtmat, family, corr, 1:15, 1:4, wgtmat,
#'       test = "dunnett", nthreads = 1)
#'
#' @export
fadjp <- function(stg2_p, wgtmat = NULL, family = NULL, corr = NULL,
                  stg1_inthyp_nr, stg2_elemhyp, stg2_wgtmat = NULL,
                  test = "dunnett", nthreads = 0) {

  m <- if (!is.null(wgtmat)) {
    ncol(wgtmat$wgtmat)
  } else if (!is.null(family)) {
    ncol(family)
  } else if (!is.null(corr)) {
    ncol(corr)
  } else {
    stop("At least one of wgtmat, family, or corr must be provided")
  }
  if (is.null(family)) {
    family <- matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family <- matrix(family, ncol = m)
  }
  if (is.null(wgtmat)) {
    wgtmat <- fDefaultWgtmat(m)
  }
  if (is.null(stg2_wgtmat)) {
    stg2_wgtmat <- fDefaultWgtmat(length(stg2_elemhyp))
  }
  if (is.null(corr)) {
    corr <- matrix(NA_real_, m, m)
    diag(corr) <- 1
    for (h in seq_len(nrow(family))) {
      idx <- which(family[h, ] != 0)
      corr[idx, idx] <- 0.5
    }
    diag(corr) <- 1
  }

  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  fadjpRcpp(stg2_p = stg2_p, wgtmat = wgtmat, family = family, corr = corr,
            stg1_inthyp_nr = stg1_inthyp_nr,  stg2_elemhyp = stg2_elemhyp,
            stg2_wgtmat = stg2_wgtmat, test = test)
}


#' @title Determine stage 1 rejections for adaptive multiple testing
#' @description Determine the elementary hypotheses rejected at stage 1 and
#'   the intersection hypotheses retained for stage 2.
#' @param m Number of elementary hypotheses.
#' @param stg1_loc_p A list returned by \code{fadjp}, containing the stage 1
#'   intersection hypotheses and local p-values.
#' @param alpha1 Stage 1 significance level.
#' @return A list containing \code{stg1_elemhyp_r_idx}, the rejected
#'   elementary hypotheses, \code{stg1_inthyp_nr_idx}, the unrejected
#'   intersection hypotheses, and \code{inthyp}, their indicator matrix.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' wgtmat <- fwgtmat(c(0.5, 0.5, 0, 0),
#'                   matrix(c(0, 0.5, 0.5, 0,
#'                            0.5, 0, 0, 0.5,
#'                            0, 1, 0, 0,
#'                            1, 0, 0, 0),
#'                          nrow = 4, byrow = TRUE))
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' stage1_local_pvalues <- fadjp(
#'   c(0.00045, 0.0952, 0.0225, 0.1104),
#'   wgtmat, family, corr,
#'   1:15, 1:4, wgtmat)
#' fPCStage1(4, stage1_local_pvalues,
#'           errorSpent(0.5, 0.025, "sfOF"))
#'
#' @export
fPCStage1 <- function(m, stg1_loc_p, alpha1) {
  fPCStage1Rcpp(m = m, stg1_loc_p = stg1_loc_p, alpha1 = alpha1)
}


#' @title Determine stage 2 rejections using p-value combination
#' @description Combine stage 1 and stage 2 local p-values and determine the
#'   elementary hypotheses rejected after stage 2.
#' @param stg1_elemhyp_r_idx Indices of the elementary hypotheses rejected at
#'   stage 1.
#' @param stg2_elemhyp_idx Indices of the elementary hypotheses tested at stage
#'   2.
#' @param stg1_loc_p A list returned by \code{fadjp} for stage 1.
#' @param stg2_loc_p A list returned by \code{fadjp} for stage 2.
#' @param alpha Overall significance level.
#' @param info_frac Information fraction for stage 1.
#' @return A list containing the combined local p-values and \code{rej_elem},
#'   an integer vector indicating the elementary hypotheses rejected after
#'   stage 2.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' alpha <- 0.025
#' alpha1 <- errorSpent(0.5, alpha, "sfOF")
#' stage1_pvalues <- c(0.00045, 0.0952, 0.0225, 0.1104)
#' stage1_loc_p <- fadjp(stage1_pvalues, wgtmat, family, corr,
#'                       1:15, 1:4, wgtmat)
#' stage1_rejections <- fPCStage1(4, stage1_loc_p, alpha1)
#'
#' adapted_graph <- updateGraph(initial_weights,
#'                              transition_matrix,
#'                              I = 1:4, j = 1)
#' stage2_weight_matrix <- fwgtmat(
#'   adapted_graph$w[adapted_graph$I],
#'   adapted_graph$G[adapted_graph$I, adapted_graph$I])
#' stage2_pvalues <- c(0.1121, 0.0112, 0.1153)
#' stage2_loc_p <- fadjp(stage2_pvalues, wgtmat, family, corr,
#'                       stage1_rejections$stg1_inthyp_nr_idx,
#'                       adapted_graph$I, stage2_weight_matrix)
#' fPCrej(stage1_rejections$stg1_elemhyp_r_idx,
#'        adapted_graph$I, stage1_loc_p,
#'        stage2_loc_p, alpha, 0.5)
#'
#' # Change weights for the elementary hypotheses.
#' stage2_weight_matrix_reweighted <- fwgtmat(
#'   w = c(0.5, 0.25, 0.25),
#'   G = adapted_graph$G[adapted_graph$I, adapted_graph$I])
#' stage2_loc_p_reweighted <- fadjp(
#'   stage2_pvalues, wgtmat, family, corr,
#'   stage1_rejections$stg1_inthyp_nr_idx, adapted_graph$I,
#'   stage2_weight_matrix_reweighted)
#' fPCrej(stage1_rejections$stg1_elemhyp_r_idx,
#'        adapted_graph$I, stage1_loc_p,
#'        stage2_loc_p_reweighted, alpha, 0.5)
#'
#' @export
fPCrej <- function(stg1_elemhyp_r_idx, stg2_elemhyp_idx, stg1_loc_p,
                   stg2_loc_p, alpha, info_frac) {
  fPCrejRcpp(stg1_elemhyp_r_idx = stg1_elemhyp_r_idx,
             stg2_elemhyp_idx = stg2_elemhyp_idx,
             stg1_loc_p = stg1_loc_p, stg2_loc_p = stg2_loc_p,
             alpha = alpha, info_frac = info_frac)
}


#' @title Compute stage-wise bounds for multiple testing
#' @description Compute stage-wise bounds for a two-stage multiple testing
#' procedure.
#' @param m Number of hypotheses
#' @param wgtmat Weight matrix for the intersection hypotheses
#' @param family Family matrix indicating which hypotheses belong to which
#'   families. The correlation is known only for hypotheses belonging to the
#'   same family. If NULL, all hypotheses are assumed to belong to the same
#'   family.
#' @param corr Correlation matrix for the test statistics
#' @param alpha Overall significance level
#' @param alpha1 Significance level for the first stage
#' @param info_frac Information fraction for the first stage
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list containing the stage 1 and stage 2 bounds for the individual
#'   hypotheses in each intersection hypothesis.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' alpha <- 0.025
#' alpha1 <- errorSpent(0.5, alpha, "sfOF")
#' fStageBound(4, wgtmat, family, corr, alpha, alpha1, 0.5,
#'             nthreads = 1)
#'
#' @export
fStageBound <- function(m, wgtmat, family = NULL, corr = NULL, alpha, alpha1,
                        info_frac, nthreads = 0) {
  if (is.null(family)) {
    family <- matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family <- matrix(family, ncol = m)
  }

  if (is.null(corr)) {
    corr <- 0.5*diag(m) + 0.5
  }

  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  x <- fStageBoundRcpp(m = m, wgtmat = wgtmat, family = family, corr = corr,
                       alpha = alpha, alpha1 = alpha1, info_frac = info_frac)
  x
}

#' @title Compute conditional error rates
#' @description Compute conditional error rates for the intersection hypotheses
#'   after the first stage of a two-stage multiple testing procedure.
#' @param m Number of hypotheses.
#' @param wgtmat Weight matrix for the intersection hypotheses.
#' @param family Family matrix indicating which hypotheses belong to which
#'   families. If NULL, all hypotheses are assumed to belong to the same family.
#' @param corr Correlation matrix for the test statistics.
#' @param info_frac Information fraction for the first stage.
#' @param stg1_bnd Stage 1 bounds for the elementary hypotheses in each
#'   intersection hypothesis.
#' @param stg2_bnd Stage 2 bounds for the elementary hypotheses in each
#'   intersection hypothesis.
#' @param stg1_p P-values for the elementary hypotheses in stage 1.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list containing the conditional error rates and the hypotheses
#'   rejected or not rejected at stage 1.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' bounds <- fStageBound(4, wgtmat, family, corr, 0.025,
#'                       errorSpent(0.5, 0.025, "sfOF"),
#'                       0.5, nthreads = 1)
#' fCER(4, wgtmat, family, corr, 0.5, bounds$stg1_bnd,
#'      bounds$stg2_bnd, c(0.00045, 0.0952, 0.0225, 0.1104),
#'      nthreads = 1)
#'
#' @export
fCER <- function(m, wgtmat, family = NULL, corr = NULL, info_frac, stg1_bnd,
                 stg2_bnd, stg1_p, nthreads = 0) {
  if (is.null(family)) {
    family <- matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family <- matrix(family, ncol = m)
  }

  if (is.null(corr)) {
    corr <- 0.5*diag(m) + 0.5
  }

  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  x <- fCERRcpp(m = m, wgtmat = wgtmat, family = family, corr = corr,
                info_frac = info_frac, stg1_bnd = stg1_bnd,
                stg2_bnd = stg2_bnd, stg1_p = stg1_p)
  x
}

#' @title Compute new stage-wise bounds after adaptation
#' @description Compute new stage 2 bounds after an adaptation of the
#'   hypotheses tested or their weights.
#' @param m Number of hypotheses.
#' @param wgtmat Weight matrix for the original intersection hypotheses.
#' @param family Family matrix indicating which hypotheses belong to which
#'   families. If NULL, all hypotheses are assumed to belong to the same family.
#' @param corr Correlation matrix for the test statistics.
#' @param stg1_p P-values for the elementary hypotheses in stage 1.
#' @param stg1_inthyp_nr_idx Indices of the intersection hypotheses not rejected
#'   in stage 1.
#' @param CER Conditional error rates for the intersection hypotheses.
#' @param stg2_elemhyp_idx Indices of the elementary hypotheses in stage 2.
#' @param stg2_wgtmat Weight matrix for the intersection hypotheses in stage 2.
#' @param info_frac_new New information fraction for stage 1 after adaptation.
#' @param nthreads The number of threads to use in simulations (0 means
#'   the default RcppParallel behavior).
#'
#' @return A list containing the adapted intersection hypotheses and their new
#'   stage 2 bounds.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1),
#'                  nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1),
#'                nrow = 4, byrow = TRUE)
#' stage1_pvalues <- c(0.00045, 0.0952, 0.0225, 0.1104)
#' bounds <- fStageBound(4, wgtmat, family, corr, 0.025,
#'                       errorSpent(0.5, 0.025, "sfOF"), 0.5)
#' conditional_error_rates <- fCER(
#'   4, wgtmat, family, corr, 0.5,
#'   bounds$stg1_bnd, bounds$stg2_bnd, stage1_pvalues)
#' stage2_weight_matrix <- fwgtmat(
#'   w = c(0.5, 0.5),
#'   G = matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
#' fNewBound(4, wgtmat, family, corr, stage1_pvalues,
#'           conditional_error_rates$stg1_inthyp_nr_idx,
#'           conditional_error_rates$CER, c(2, 4),
#'           stage2_weight_matrix, 0.4, nthreads = 1)
#'
#' @export
fNewBound <- function(m, wgtmat, family = NULL, corr = NULL, stg1_p,
                      stg1_inthyp_nr_idx, CER, stg2_elemhyp_idx,
                      stg2_wgtmat, info_frac_new, nthreads = 0) {
  if (is.null(family)) {
    family <- matrix(1, 1, m)
  } else if (!is.matrix(family)) {
    family <- matrix(family, ncol = m)
  }

  if (is.null(corr)) {
    corr <- 0.5*diag(m) + 0.5
  }

  if (nthreads > 0) {
    n_physical_cores <- parallel::detectCores(logical = FALSE)
    RcppParallel::setThreadOptions(min(nthreads, n_physical_cores))
  }

  x <- fNewBoundRcpp(m = m, wgtmat = wgtmat, family = family, corr = corr,
                     stg1_p = stg1_p, stg1_inthyp_nr_idx = stg1_inthyp_nr_idx,
                     CER = CER, stg2_elemhyp_idx = stg2_elemhyp_idx,
                     stg2_wgtmat = stg2_wgtmat, info_frac_new = info_frac_new)
  x
}

#' @title Determine stage 2 rejections from conditional error rates
#' @description Determine the elementary hypotheses rejected after stage 2
#'   using adapted stage 2 bounds and cumulative p-values.
#' @param stg1_elemhyp_r_idx Indices of the elementary hypotheses rejected
#'   in stage 1.
#' @param stg2_elemhyp_idx Indices of the elementary hypotheses in stage 2.
#' @param stg2_inthyp_idx Indices of the intersection hypotheses in stage 2.
#' @param stg2_bnd_new New stage 2 bounds for the elementary hypotheses in each
#'   intersection hypothesis.
#' @param cum_p Cumulative p-values for the elementary hypotheses in stage 2.
#' @return A logical vector indicating whether each individual hypothesis is
#'   rejected after stage 2.
#'
#' @references
#' Cyrus Mehta, Ajoy Mukhopadhyay, and Martin Posch. Graph Based, Adaptive,
#' Multiarm, Multiple Endpoint, Two-Stage Designs. Statistics in Medicine.
#' 2025.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @examples
#' initial_weights <- c(0.5, 0.5, 0, 0)
#' transition_matrix <- matrix(c(0, 0.5, 0.5, 0,
#'                               0.5, 0, 0, 0.5,
#'                               0, 1, 0, 0,
#'                               1, 0, 0, 0),
#'                             nrow = 4, byrow = TRUE)
#' wgtmat <- fwgtmat(initial_weights, transition_matrix)
#' family <- matrix(c(1, 1, 0, 0,
#'                    0, 0, 1, 1), nrow = 2, byrow = TRUE)
#' corr <- matrix(c(1, 0.5, NA, NA,
#'                  0.5, 1, NA, NA,
#'                  NA, NA, 1, 0.5,
#'                  NA, NA, 0.5, 1), nrow = 4, byrow = TRUE)
#' stage1_pvalues <- c(0.00045, 0.0952, 0.0225, 0.1104)
#' bounds <- fStageBound(4, wgtmat, family, corr, 0.025,
#'                       errorSpent(0.5, 0.025, "sfOF"), 0.5)
#' conditional_error_rates <- fCER(
#'   4, wgtmat, family, corr, 0.5,
#'   bounds$stg1_bnd, bounds$stg2_bnd, stage1_pvalues)
#' stage2_weight_matrix <- fwgtmat(
#'   w = c(0.5, 0.5),
#'   G = matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE))
#' adapted_bounds <- fNewBound(
#'   4, wgtmat, family, corr, stage1_pvalues,
#'   conditional_error_rates$stg1_inthyp_nr_idx,
#'   conditional_error_rates$CER,
#'   c(2, 4), stage2_weight_matrix, 0.4)
#' stage2_cumulative_pvalues <-
#'   1 - pnorm(sqrt(0.4) * qnorm(1 - c(0.0952, 0.1104)) +
#'               sqrt(1 - 0.4) * qnorm(1 - c(0.0299, 0.0586)))
#' fCERrej(conditional_error_rates$stg1_elemhyp_r_idx, c(2, 4),
#'         adapted_bounds$inthyp, adapted_bounds$stg2_bnd_new,
#'         stage2_cumulative_pvalues)
#'
#' @export
fCERrej <- function(stg1_elemhyp_r_idx, stg2_elemhyp_idx, stg2_inthyp_idx,
                    stg2_bnd_new, cum_p) {
  fCERrejRcpp(stg1_elemhyp_r_idx, stg2_elemhyp_idx, stg2_inthyp_idx,
              stg2_bnd_new, cum_p)
}

