testthat::test_that("CI endpoints are ordered and contain point estimates where expected", {
  cp <- ClopperPearsonCI(n = 20, y = 3)
  testthat::expect_equal(round(cp$phat, 6), 0.15)
  testthat::expect_true(cp$lower <= cp$phat && cp$phat <= cp$upper)
  testthat::expect_equal(round(cp$lower, 6), 0.032071)
  testthat::expect_equal(round(cp$upper, 6), 0.378927)

  mn_or <- mnOddsRatioCI(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0))
  mn_rd <- mnRiskDiffCI(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0))
  mn_rr <- mnRiskRatioCI(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0))
  mn_rate_d <- mnRateDiffCI(t1 = c(10, 10), y1 = c(4, 3), t2 = c(20, 10), y2 = c(2, 0))
  mn_rate_r <- mnRateRatioCI(t1 = c(10, 10), y1 = c(4, 3), t2 = c(20, 10), y2 = c(2, 0))

  for (ci in list(mn_or, mn_rd, mn_rr, mn_rate_d, mn_rate_r)) {
    testthat::expect_true(ci$lower <= ci$estimate && ci$estimate <= ci$upper)
  }

  testthat::expect_equal(mn_or$estimate, 12.95689, tolerance = 1e-5)
  testthat::expect_equal(round(mn_or$lower, 6), 1.682969)
  testthat::expect_equal(round(mn_or$upper, 4), 158.8478)

  testthat::expect_equal(round(mn_rd$estimate, 6), 0.3)
  testthat::expect_equal(round(mn_rd$lower, 6), 0.082536)
  testthat::expect_equal(round(mn_rd$upper, 5), 0.53206)

  testthat::expect_equal(round(mn_rr$estimate, 2), 6.25)
  testthat::expect_equal(round(mn_rr$lower, 6), 1.475588)
  testthat::expect_equal(round(mn_rr$upper, 4), 27.6909)

  testthat::expect_equal(round(mn_rate_d$estimate, 6), 0.3)
  testthat::expect_equal(round(mn_rate_d$lower, 6), 0.065117)
  testthat::expect_equal(round(mn_rate_d$upper, 4), 0.6871)

  testthat::expect_equal(round(mn_rate_r$estimate, 2), 6.25)
  testthat::expect_equal(round(mn_rate_r$lower, 6), 1.339243)
  testthat::expect_equal(mn_rate_r$upper, 29.16759, tolerance = 1e-4)

  rd_exact <- riskDiffExactCI(n1 = 30, y1 = 2, n2 = 30, y2 = 1, cilevel = 0.95)
  rr_exact <- riskRatioExactCI(n1 = 30, y1 = 2, n2 = 30, y2 = 1, cilevel = 0.95)
  testthat::expect_true(rd_exact$lower <= rd_exact$estimate && rd_exact$estimate <= rd_exact$upper)
  testthat::expect_true(rr_exact$lower <= rr_exact$estimate && rr_exact$estimate <= rr_exact$upper)
  testthat::expect_equal(round(rd_exact$estimate, 6), 0.033333)
  testthat::expect_equal(round(rd_exact$lower, 6), -0.111153)
  testthat::expect_equal(round(rd_exact$upper, 6), 0.200071)
  testthat::expect_equal(round(rr_exact$estimate, 0), 2)
  testthat::expect_equal(round(rr_exact$lower, 6), 0.185541)
  testthat::expect_equal(round(rr_exact$upper, 5), 54.22365)

  ro <- remlOddsRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, oddsRatioH0 = 1.25)
  rd <- remlRiskDiff(n1 = 10, y1 = 4, n2 = 20, y2 = 0, riskDiffH0 = 0.1)
  rr <- remlRiskRatio(n1 = 10, y1 = 4, n2 = 20, y2 = 2, riskRatioH0 = 1.2)
  rtd <- remlRateDiff(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateDiffH0 = 0.1)
  rtr <- remlRateRatio(t1 = 10, y1 = 4, t2 = 20, y2 = 2, rateRatioH0 = 1.1)

  testthat::expect_equal(round((ro[1] / (1 - ro[1])) / (ro[2] / (1 - ro[2])), 6), 1.25)
  testthat::expect_equal(round(rd[1] - rd[2], 6), 0.1)
  testthat::expect_equal(round(rr[1] / rr[2], 6), 1.2)
  testthat::expect_equal(round(rtd[1] - rtd[2], 6), 0.1)
  testthat::expect_equal(round(rtr[1] / rtr[2], 6), 1.1)

  testthat::expect_equal(round(zstatOddsRatio(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0), oddsRatioH0 = 1), 6), 2.64148)
  testthat::expect_equal(round(zstatRiskDiff(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0), riskDiffH0 = 0), 6), 2.627423)
  testthat::expect_equal(round(zstatRiskRatio(n1 = c(10, 10), y1 = c(4, 3), n2 = c(20, 10), y2 = c(2, 0), riskRatioH0 = 1), 6), 2.627423)
  testthat::expect_equal(round(zstatRateDiff(t1 = c(10, 10), y1 = c(4, 3), t2 = c(20, 10), y2 = c(2, 0), rateDiffH0 = 0), 6), 2.424871)
  testthat::expect_equal(round(zstatRateRatio(t1 = c(10, 10), y1 = c(4, 3), t2 = c(20, 10), y2 = c(2, 0), rateRatioH0 = 1), 6), 2.424871)
})


testthat::test_that("p-values agree with CI inversion on selected points", {
  rd_ci <- riskDiffExactCI(n1 = 68, y1 = 2, n2 = 65, y2 = 1, cilevel = 0.95)
  rr_ci <- riskRatioExactCI(n1 = 68, y1 = 2, n2 = 65, y2 = 1, cilevel = 0.95)

  rd_null <- riskDiffExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskDiffH0 = 0, directionUpper = TRUE)
  rr_null <- riskRatioExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskRatioH0 = 1, directionUpper = TRUE)

  rd_at_lower <- riskDiffExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskDiffH0 = rd_ci[["lower"]][1], directionUpper = TRUE)
  rr_at_lower <- riskRatioExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskRatioH0 = rr_ci[["lower"]][1], directionUpper = TRUE)
  rd_at_upper <- riskDiffExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskDiffH0 = rd_ci[["upper"]][1], directionUpper = TRUE)
  rr_at_upper <- riskRatioExactPValue(n1 = 68, y1 = 2, n2 = 65, y2 = 1, riskRatioH0 = rr_ci[["upper"]][1], directionUpper = TRUE)

  testthat::expect_equal(round(rd_null$pvalue, 6), 0.354387)
  testthat::expect_equal(round(rr_null$pvalue, 6), 0.354387)

  # For one-sided tests, lower CI bound corresponds to p ≈ alpha/2.
  testthat::expect_equal(round(rd_at_lower$pvalue, 3), 0.025)
  testthat::expect_equal(round(rr_at_lower$pvalue, 3), 0.025)

  # Upper CI bound gives near-one upper-tail p-value.
  testthat::expect_gt(rd_at_upper$pvalue, 0.95)
  testthat::expect_gt(rr_at_upper$pvalue, 0.95)
})


testthat::test_that("zero-cell and all-event edge cases are valid", {
  cp_zero <- ClopperPearsonCI(20, 0)
  cp_all <- ClopperPearsonCI(20, 20)
  testthat::expect_equal(cp_zero$lower, 0)
  testthat::expect_equal(cp_all$upper, 1)
  testthat::expect_true(cp_zero$upper < 0.2)
  testthat::expect_true(cp_all$lower > 0.8)

  mn_rr_zero <- mnRiskRatioCI(n1 = 20, y1 = 0, n2 = 20, y2 = 2)
  mn_or_zero <- mnOddsRatioCI(n1 = 20, y1 = 0, n2 = 20, y2 = 2)
  mn_rd_all <- mnRiskDiffCI(n1 = 20, y1 = 20, n2 = 20, y2 = 20)
  testthat::expect_equal(mn_rr_zero$estimate, 0)
  testthat::expect_equal(mn_or_zero$estimate, 0)
  testthat::expect_true(mn_rr_zero$upper > 1)
  testthat::expect_true(mn_or_zero$upper > 1)
  testthat::expect_true(mn_rd_all$lower <= 0 && 0 <= mn_rd_all$upper)

  rd_ci_zero <- riskDiffExactCI(n1 = 20, y1 = 0, n2 = 20, y2 = 2)
  rr_ci_all <- riskRatioExactCI(n1 = 20, y1 = 20, n2 = 20, y2 = 20)
  testthat::expect_true(rd_ci_zero$lower <= rd_ci_zero$estimate && rd_ci_zero$estimate <= rd_ci_zero$upper)
  testthat::expect_true(rr_ci_all$lower <= rr_ci_all$estimate && rr_ci_all$estimate <= rr_ci_all$upper)

  rd_p_zero <- riskDiffExactPValue(n1 = 20, y1 = 0, n2 = 20, y2 = 2, riskDiffH0 = 0, directionUpper = TRUE)
  rr_p_all <- riskRatioExactPValue(n1 = 20, y1 = 20, n2 = 20, y2 = 20, riskRatioH0 = 1, directionUpper = TRUE)
  testthat::expect_true(rd_p_zero$pvalue >= 0 && rd_p_zero$pvalue <= 1)
  testthat::expect_true(rr_p_all$pvalue >= 0 && rr_p_all$pvalue <= 1)

  testthat::expect_true(is.finite(zstatRiskRatio(n1 = 20, y1 = 0, n2 = 20, y2 = 2, riskRatioH0 = 1)))
  testthat::expect_true(is.finite(zstatOddsRatio(n1 = 20, y1 = 0, n2 = 20, y2 = 2, oddsRatioH0 = 1)))
})