adaptive_two_stage_setup <- function() {
  initial_weights <- c(0.5, 0.5, 0, 0)
  transition_matrix <- matrix(
    c(0, 0.5, 0.5, 0,
      0.5, 0, 0, 0.5,
      0, 1, 0, 0,
      1, 0, 0, 0),
    nrow = 4,
    byrow = TRUE
  )
  wgtmat <- fwgtmat(initial_weights, transition_matrix)
  family <- matrix(
    c(1, 1, 0, 0,
      0, 0, 1, 1),
    nrow = 2,
    byrow = TRUE
  )
  corr <- matrix(
    c(1, 0.5, NA, NA,
      0.5, 1, NA, NA,
      NA, NA, 1, 0.5,
      NA, NA, 0.5, 1),
    nrow = 4,
    byrow = TRUE
  )
  list(
    initial_weights = initial_weights,
    transition_matrix = transition_matrix,
    wgtmat = wgtmat,
    family = family,
    corr = corr
  )
}


testthat::test_that("fadjp returns stage-wise local p-values", {
  setup <- adaptive_two_stage_setup()
  local_pvalues <- lrstat:::fadjp(
    c(0.00045, 0.0952, 0.0225, 0.1104), setup$wgtmat, setup$family,
    setup$corr, 1:15, 1:4, setup$wgtmat
  )

  default_corr_pvalues <- lrstat:::fadjp(
    c(0.00045, 0.0952, 0.0225, 0.1104), wgtmat = setup$wgtmat,
    family = setup$family, stg1_inthyp_nr = 1:15, stg2_elemhyp = 1:4,
    stg2_wgtmat = setup$wgtmat
  )

  testthat::expect_named(local_pvalues, c("inthyp_idx", "inthyp", "pinter"))
  testthat::expect_equal(local_pvalues$inthyp_idx, seq_len(15))
  testthat::expect_equal(dim(local_pvalues$inthyp), c(15L, 4L))
  testthat::expect_length(local_pvalues$pinter, 15L)
  testthat::expect_true(all(local_pvalues$pinter >= 0 & local_pvalues$pinter <= 1))
  testthat::expect_equal(default_corr_pvalues$pinter, local_pvalues$pinter)
})


testthat::test_that("fPCStage1 identifies stage 1 rejections and retained intersections", {
  setup <- adaptive_two_stage_setup()
  local_pvalues <- lrstat:::fadjp(
    c(0.00045, 0.0952, 0.0225, 0.1104), setup$wgtmat, setup$family,
    setup$corr, 1:15, 1:4, setup$wgtmat
  )
  stage1 <- lrstat:::fPCStage1(4, local_pvalues, errorSpent(0.5, 0.025, "sfOF"))

  testthat::expect_named(
    stage1,
    c("stg1_elemhyp_r_idx", "stg1_inthyp_nr_idx", "inthyp")
  )
  testthat::expect_true(all(stage1$stg1_elemhyp_r_idx %in% 1:4))
  testthat::expect_true(all(stage1$stg1_inthyp_nr_idx %in% 1:15))
  testthat::expect_equal(ncol(stage1$inthyp), 4L)
  testthat::expect_equal(
    nrow(stage1$inthyp),
    length(stage1$stg1_inthyp_nr_idx)
  )
})


testthat::test_that("fStageBound returns valid stage-wise boundary matrices", {
  setup <- adaptive_two_stage_setup()
  stage_bounds <- lrstat:::fStageBound(
    4, setup$wgtmat, setup$family, setup$corr,
    alpha = 0.025, alpha1 = errorSpent(0.5, 0.025, "sfOF"),
    info_frac = 0.5
  )

  testthat::expect_named(
    stage_bounds,
    c("inthyp", "stg1_coef", "stg2_coef", "stg1_bnd", "stg2_bnd")
  )
  testthat::expect_equal(dim(stage_bounds$inthyp), c(15L, 4L))
  testthat::expect_equal(dim(stage_bounds$stg1_bnd), c(15L, 4L))
  testthat::expect_equal(dim(stage_bounds$stg2_bnd), c(15L, 4L))
  testthat::expect_true(all(stage_bounds$stg1_bnd >= 0 & stage_bounds$stg1_bnd <= 1))
  testthat::expect_true(all(stage_bounds$stg2_bnd >= 0 & stage_bounds$stg2_bnd <= 1))
  testthat::expect_error(
    lrstat:::fStageBound(4, setup$wgtmat, setup$family, setup$corr,
                alpha = 0.025, alpha1 = 0.025, info_frac = 0.5),
    regexp = "alpha1 must be in"
  )
})


testthat::test_that("fCER computes conditional error rates from stage bounds", {
  setup <- adaptive_two_stage_setup()
  stage_bounds <- lrstat:::fStageBound(
    4, setup$wgtmat, setup$family, setup$corr,
    alpha = 0.025, alpha1 = errorSpent(0.5, 0.025, "sfOF"),
    info_frac = 0.5
  )
  conditional_error <- lrstat:::fCER(
    4, setup$wgtmat, setup$family, setup$corr, 0.5,
    stage_bounds$stg1_bnd, stage_bounds$stg2_bnd,
    c(0.00045, 0.0952, 0.0225, 0.1104)
  )

  testthat::expect_named(
    conditional_error,
    c("stg1_elemhyp_r_idx", "stg1_inthyp_nr_idx", "inthyp", "CER")
  )
  testthat::expect_equal(ncol(conditional_error$inthyp), 4L)
  testthat::expect_equal(nrow(conditional_error$inthyp), length(conditional_error$CER))
  testthat::expect_true(all(conditional_error$CER >= 0 & conditional_error$CER <= 1))
})


testthat::test_that("fPCrej combines stage p-values and preserves rejection decisions", {
  setup <- adaptive_two_stage_setup()
  stage1_local_pvalues <- lrstat:::fadjp(
    c(0.00045, 0.0952, 0.0225, 0.1104), setup$wgtmat, setup$family,
    setup$corr, 1:15, 1:4, setup$wgtmat
  )
  stage1 <- lrstat:::fPCStage1(4, stage1_local_pvalues, errorSpent(0.5, 0.025, "sfOF"))
  adapted_graph <- updateGraph(
    setup$initial_weights, setup$transition_matrix, I = 1:4, j = 1
  )
  stage2_weight_matrix <- fwgtmat(
    adapted_graph$w[adapted_graph$I],
    adapted_graph$G[adapted_graph$I, adapted_graph$I]
  )
  stage2_local_pvalues <- lrstat:::fadjp(
    c(0.1121, 0.0112, 0.1153), setup$wgtmat, setup$family, setup$corr,
    stage1$stg1_inthyp_nr_idx, adapted_graph$I, stage2_weight_matrix
  )
  combined <- lrstat:::fPCrej(
    stage1$stg1_elemhyp_r_idx, adapted_graph$I,
    stage1_local_pvalues, stage2_local_pvalues, 0.025, 0.5
  )

  testthat::expect_named(
    combined,
    c("stg1_inthyp_nr_idx", "stg2_elemhyp_idx", "inthyp",
      "stg1_pinter", "stg2_pinter", "cum_pinter", "rej_elem")
  )
  testthat::expect_equal(combined$stg2_elemhyp_idx, adapted_graph$I)
  testthat::expect_equal(length(combined$rej_elem), 4L)
  testthat::expect_true(all(combined$rej_elem %in% 0:1))
  testthat::expect_true(all(combined$cum_pinter >= 0 & combined$cum_pinter <= 1))
})


testthat::test_that("fNewBound and fCERrej complete the conditional-error workflow", {
  setup <- adaptive_two_stage_setup()
  stage_bounds <- lrstat:::fStageBound(
    4, setup$wgtmat, setup$family, setup$corr,
    alpha = 0.025, alpha1 = errorSpent(0.5, 0.025, "sfOF"),
    info_frac = 0.5
  )
  stage1_pvalues <- c(0.00045, 0.0952, 0.0225, 0.1104)
  conditional_error <- lrstat:::fCER(
    4, setup$wgtmat, setup$family, setup$corr, 0.5,
    stage_bounds$stg1_bnd, stage_bounds$stg2_bnd, stage1_pvalues
  )
  stage2_weight_matrix <- fwgtmat(
    c(0.5, 0.5), matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE)
  )
  adapted_bounds <- lrstat:::fNewBound(
    4, setup$wgtmat, setup$family, setup$corr, stage1_pvalues,
    conditional_error$stg1_inthyp_nr_idx, conditional_error$CER,
    c(2, 4), stage2_weight_matrix, 0.4
  )
  cumulative_pvalues <- 1 - pnorm(
    sqrt(0.4) * qnorm(1 - c(0.0952, 0.1104)) +
      sqrt(1 - 0.4) * qnorm(1 - c(0.0299, 0.0586))
  )
  rejection <- lrstat:::fCERrej(
    conditional_error$stg1_elemhyp_r_idx, c(2, 4),
    adapted_bounds$inthyp, adapted_bounds$stg2_bnd_new,
    cumulative_pvalues
  )

  testthat::expect_named(adapted_bounds, c("inthyp", "stg2_coef_new", "stg2_bnd_new"))
  testthat::expect_equal(ncol(adapted_bounds$inthyp), 4L)
  testthat::expect_equal(dim(adapted_bounds$stg2_bnd_new), dim(adapted_bounds$inthyp))
  testthat::expect_true(all(adapted_bounds$stg2_bnd_new >= 0 & adapted_bounds$stg2_bnd_new <= 1))
  testthat::expect_type(rejection, "logical")
  testthat::expect_length(rejection, 4L)
})
