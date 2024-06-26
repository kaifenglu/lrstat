testthat::test_that("accrual duration with specified critical values, unweighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   criticalValues = c(2.250, 2.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = NA,
                   followupTime = 18, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 0,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$accrualDuration, 2),
                         23.56)
})


testthat::test_that("accrual duration with specified critical values, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   criticalValues = c(2.250, 2.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = NA,
                   followupTime = 18, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$accrualDuration, 2),
                         16.74)
})


testthat::test_that("accrual duration with alpha-spending, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   alpha = 0.025, typeAlphaSpending = "user",
                   userAlphaSpending = c(0.01, 0.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = NA,
                   followupTime = 18, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$accrualDuration, 2),
                         16.67)
})


testthat::test_that("accrual duration with beta-spending, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   alpha = 0.025, typeAlphaSpending = "user",
                   userAlphaSpending = c(0.01, 0.025),
                   typeBetaSpending = "user",
                   userBetaSpending = c(0.1, 0.2),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = NA,
                   followupTime = 18, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$accrualDuration, 2),
                         16.94)
})


testthat::test_that("accrual duration with user-specified alpha and beta spending", {
  l = lrsamplesize(beta=0.1, kMax=3,
                   informationRates=c(0.5, 0.75, 1),
                   alpha=0.025, typeAlphaSpending = "user",
                   userAlphaSpending = c(0, 0.00965, 0.025),
                   typeBetaSpending = "user",
                   userBetaSpending = c(0.02, 0.02, 0.1),
                   allocationRatioPlanned = 3,
                   accrualIntensity = 5,
                   lambda1 = 0.3*0.95/12,
                   lambda2 = 0.95/12,
                   gamma1 = -log(1-0.1)/24,
                   gamma2 = -log(1-0.1)/24,
                   accrualDuration = NA,
                   followupTime = 6.5, fixedFollowup = TRUE,
                   rho1 = 0, rho2 = 0,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$accrualDuration, 2),
                         26.17)
})


testthat::test_that("follow-up duration with specified critical values, unweighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   criticalValues = c(2.250, 2.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = 22,
                   followupTime = NA, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 0,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$followupTime, 2),
                         21.43)
})


testthat::test_that("follow-up duration with specified critical values, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   criticalValues = c(2.250, 2.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = 22,
                   followupTime = NA, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$followupTime, 2),
                         9.19)
})


testthat::test_that("follow-up duration with alpha-spending, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   alpha = 0.025, typeAlphaSpending = "user",
                   userAlphaSpending = c(0.01, 0.025),
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = 22,
                   followupTime = NA, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$followupTime, 2),
                         9.06)
})


testthat::test_that("follow-up duration with alpha and beta spending, weighted", {
  l = lrsamplesize(beta = 0.2, kMax = 2,
                   informationRates = c(0.8, 1),
                   alpha = 0.025, typeAlphaSpending = "sfKD",
                   parameterAlphaSpending = 0.5,
                   typeBetaSpending = "sfHSD",
                   parameterBetaSpending = -2,
                   allocationRatioPlanned = 1,
                   accrualTime = seq(0, 8),
                   accrualIntensity = 26/9*seq(1, 9),
                   piecewiseSurvivalTime = c(0, 6),
                   stratumFraction = c(0.2, 0.8),
                   lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                   gamma1 = -log(1-0.05)/12,
                   gamma2 = -log(1-0.05)/12,
                   accrualDuration = 22,
                   followupTime = NA, fixedFollowup = FALSE,
                   rho1 = 0, rho2 = 1,
                   numSubintervals = 300,
                   rounding = 0)
  testthat::expect_equal(round(l$resultsUnderH1$overallResults$followupTime, 2),
                         10.94)
})

