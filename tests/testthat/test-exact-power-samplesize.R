testthat::test_that("power and sample size inversion consistency checks", {
  fisher_power <- powerFisherExact(n = 136, pi1 = 0.25, pi2 = 0.05, alpha = 0.05)
  fisher_ss <- samplesizeFisherExact(beta = 0.1, pi1 = 0.25, pi2 = 0.05, alpha = 0.05)
  testthat::expect_equal(round(fisher_power$power, 6), 0.897983)
  testthat::expect_equal(fisher_ss$n, 137)
  testthat::expect_equal(round(fisher_ss$power, 6), 0.906564)
  testthat::expect_equal(round(powerFisherExact(n = fisher_ss$n, pi1 = 0.25, pi2 = 0.05, alpha = 0.05)$power, 6), 0.906564)

  one_prop_power <- powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.05)
  one_prop_ss <- samplesizeOnePropExact(beta = 0.2, piH0 = 0.15, pi = 0.25, alpha = 0.025)
  testthat::expect_equal(round(one_prop_power$power, 6), 0.809674)
  testthat::expect_equal(one_prop_ss$n, 136)
  testthat::expect_equal(round(one_prop_ss$power, 6), 0.812676)
  testthat::expect_equal(round(powerOnePropExact(n = one_prop_ss$n, piH0 = 0.15, pi = 0.25, alpha = 0.025)$power, 6), 0.812676)

  one_rate_power <- powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012, D = 0.5, alpha = 0.025)
  one_rate_ss <- samplesizeOneRateExact(beta = 0.2, lambdaH0 = 0.2, lambda = 0.3, D = 1, alpha = 0.05)
  testthat::expect_equal(round(one_rate_power$power, 6), 0.900210)
  testthat::expect_equal(one_rate_ss$n, 162)
  testthat::expect_equal(round(one_rate_ss$power, 6), 0.807805)
  testthat::expect_equal(round(powerOneRateExact(n = one_rate_ss$n, lambdaH0 = 0.2, lambda = 0.3, D = 1, alpha = 0.05)$power, 6), 0.807805)

  risk_diff_power <- powerRiskDiffExact(n = 68, riskDiffH0 = 0, pi1 = 0.6, pi2 = 0.25, alpha = 0.025)
  risk_diff_ss <- samplesizeRiskDiffExact(beta = 0.2, riskDiffH0 = -0.3, pi1 = 0.8, pi2 = 0.8, alpha = 0.025)
  testthat::expect_equal(round(risk_diff_power$power, 6), 0.826318)
  testthat::expect_equal(risk_diff_ss$n, 60)
  testthat::expect_equal(round(risk_diff_ss$power, 6), 0.811130)
  testthat::expect_equal(round(powerRiskDiffExact(n = risk_diff_ss$n, riskDiffH0 = -0.3, pi1 = 0.8, pi2 = 0.8, alpha = 0.025)$power, 6), 0.811130)

  risk_diff_equiv_power <- powerRiskDiffExactEquiv(
    n = 200,
    riskDiffLower = -0.2,
    riskDiffUpper = 0.2,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05
  )
  risk_diff_equiv_ss <- samplesizeRiskDiffExactEquiv(
    beta = 0.2,
    riskDiffLower = -0.3,
    riskDiffUpper = 0.3,
    pi1 = 0.9,
    pi2 = 0.9,
    alpha = 0.05
  )
  testthat::expect_equal(round(risk_diff_equiv_power$power, 6), 0.914686)
  testthat::expect_equal(risk_diff_equiv_ss$n, 47)
  testthat::expect_equal(round(risk_diff_equiv_ss$power, 6), 0.800266)
  testthat::expect_equal(round(powerRiskDiffExactEquiv(n = risk_diff_equiv_ss$n, riskDiffLower = -0.3, riskDiffUpper = 0.3, pi1 = 0.9, pi2 = 0.9, alpha = 0.05)$power, 6), 0.800266)

  risk_ratio_power <- powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.025)
  risk_ratio_ss <- samplesizeRiskRatioExact(beta = 0.2, riskRatioH0 = 0.8, pi1 = 0.95, pi2 = 0.95, alpha = 0.05)
  testthat::expect_equal(round(risk_ratio_power$power, 6), 0.826318)
  testthat::expect_equal(risk_ratio_ss$n, 51)
  testthat::expect_equal(round(risk_ratio_ss$power, 6), 0.809295)
  testthat::expect_equal(round(powerRiskRatioExact(n = risk_ratio_ss$n, riskRatioH0 = 0.8, pi1 = 0.95, pi2 = 0.95, alpha = 0.05)$power, 6), 0.809295)

  risk_ratio_equiv_power <- powerRiskRatioExactEquiv(
    n = 200,
    riskRatioLower = 0.8,
    riskRatioUpper = 1.25,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05
  )
  risk_ratio_equiv_ss <- samplesizeRiskRatioExactEquiv(
    beta = 0.2,
    riskRatioLower = 0.7,
    riskRatioUpper = 1 / 0.7,
    pi1 = 0.95,
    pi2 = 0.95,
    alpha = 0.05
  )
  testthat::expect_equal(round(risk_ratio_equiv_power$power, 6), 0.751415)
  testthat::expect_equal(risk_ratio_equiv_ss$n, 37)
  testthat::expect_equal(round(risk_ratio_equiv_ss$power, 6), 0.827160)
  testthat::expect_equal(round(powerRiskRatioExactEquiv(n = risk_ratio_equiv_ss$n, riskRatioLower = 0.7, riskRatioUpper = 1 / 0.7, pi1 = 0.95, pi2 = 0.95, alpha = 0.05)$power, 6), 0.827160)
})


testthat::test_that("one-sided and two-sided option checks behave as expected", {
  fisher <- powerFisherExact(n = 136, pi1 = 0.25, pi2 = 0.05, alpha = 0.05)
  testthat::expect_equal(fisher$alpha, 0.05)

  one_prop_025 <- powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.025)
  one_prop_05 <- powerOnePropExact(n = 110, piH0 = 0.15, pi = 0.25, alpha = 0.05)
  testthat::expect_equal(one_prop_025$alpha, 0.025)
  testthat::expect_equal(one_prop_05$alpha, 0.05)
  testthat::expect_gt(one_prop_05$power, one_prop_025$power)

  one_rate_025 <- powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012, D = 0.5, alpha = 0.025)
  one_rate_05 <- powerOneRateExact(n = 525, lambdaH0 = 0.049, lambda = 0.012, D = 0.5, alpha = 0.05)
  testthat::expect_equal(one_rate_025$alpha, 0.025)
  testthat::expect_equal(one_rate_05$alpha, 0.05)
  testthat::expect_gt(one_rate_05$power, one_rate_025$power)

  risk_diff_025 <- powerRiskDiffExact(n = 68, riskDiffH0 = 0, pi1 = 0.6, pi2 = 0.25, alpha = 0.025, calculateAttainedAlpha = TRUE)
  risk_diff_05 <- powerRiskDiffExact(n = 68, riskDiffH0 = 0, pi1 = 0.6, pi2 = 0.25, alpha = 0.05, calculateAttainedAlpha = TRUE)
  risk_diff_no_attained <- powerRiskDiffExact(n = 68, riskDiffH0 = 0, pi1 = 0.6, pi2 = 0.25, alpha = 0.025, calculateAttainedAlpha = FALSE)
  testthat::expect_equal(risk_diff_025$alpha, 0.025)
  testthat::expect_equal(risk_diff_05$alpha, 0.05)
  testthat::expect_gt(risk_diff_05$power, risk_diff_025$power)
  testthat::expect_true("attainedAlpha" %in% names(risk_diff_025))
  testthat::expect_false("attainedAlpha" %in% names(risk_diff_no_attained))

  risk_ratio_025 <- powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.025, calculateAttainedAlpha = TRUE)
  risk_ratio_05 <- powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.05, calculateAttainedAlpha = TRUE)
  risk_ratio_no_attained <- powerRiskRatioExact(n = 68, pi1 = 0.6, pi2 = 0.25, alpha = 0.025, calculateAttainedAlpha = FALSE)
  testthat::expect_equal(risk_ratio_025$alpha, 0.025)
  testthat::expect_equal(risk_ratio_05$alpha, 0.05)
  testthat::expect_gt(risk_ratio_05$power, risk_ratio_025$power)
  testthat::expect_true("attainedAlpha" %in% names(risk_ratio_025))
  testthat::expect_false("attainedAlpha" %in% names(risk_ratio_no_attained))

  risk_diff_equiv_no_attained <- powerRiskDiffExactEquiv(
    n = 200,
    riskDiffLower = -0.2,
    riskDiffUpper = 0.2,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05,
    calculateAttainedAlpha = FALSE
  )
  risk_ratio_equiv_no_attained <- powerRiskRatioExactEquiv(
    n = 200,
    riskRatioLower = 0.8,
    riskRatioUpper = 1.25,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05,
    calculateAttainedAlpha = FALSE
  )
  testthat::expect_false("attainedAlphaH10" %in% names(risk_diff_equiv_no_attained))
  testthat::expect_false("attainedAlphaH10" %in% names(risk_ratio_equiv_no_attained))
})


testthat::test_that("equivalence exact methods validate margin constraints", {
  power_risk_diff_bad <- powerRiskDiffExactEquiv(
    n = 100,
    riskDiffLower = 0.2,
    riskDiffUpper = -0.2,
    pi1 = 0.8,
    pi2 = 0.8,
    alpha = 0.05
  )
  testthat::expect_equal(power_risk_diff_bad$power, 0)
  testthat::expect_equal(power_risk_diff_bad$attainedAlphaH10, 0)
  testthat::expect_equal(power_risk_diff_bad$attainedAlphaH20, 0)

  testthat::expect_error(
    samplesizeRiskDiffExactEquiv(
      beta = 0.2,
      riskDiffLower = 0.2,
      riskDiffUpper = -0.2,
      pi1 = 0.9,
      pi2 = 0.9,
      alpha = 0.05
    )
  )

  power_risk_ratio_bad <- powerRiskRatioExactEquiv(
    n = 100,
    riskRatioLower = 1.2,
    riskRatioUpper = 1.1,
    pi1 = 0.8,
    pi2 = 0.8,
    alpha = 0.05
  )
  testthat::expect_equal(power_risk_ratio_bad$power, 0)
  testthat::expect_equal(power_risk_ratio_bad$attainedAlphaH10, 0)
  testthat::expect_equal(power_risk_ratio_bad$attainedAlphaH20, 0)

  testthat::expect_error(
    samplesizeRiskRatioExactEquiv(
      beta = 0.2,
      riskRatioLower = 1.2,
      riskRatioUpper = 1.1,
      pi1 = 0.95,
      pi2 = 0.95,
      alpha = 0.05
    )
  )
})