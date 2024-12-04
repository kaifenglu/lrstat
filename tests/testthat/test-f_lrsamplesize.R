testthat::test_that(
  "accrual duration given power and follow-up time", {
    l = lrsamplesize(beta = 0.2, kMax = 2,
                     informationRates = c(0.8, 1),
                     alpha = 0.025, typeAlphaSpending = "sfOF",
                     accrualTime = seq(0, 8),
                     accrualIntensity = 26/9*seq(1, 9),
                     piecewiseSurvivalTime = c(0, 6),
                     stratumFraction = c(0.2, 0.8),
                     lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                     lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                     gamma1 = -log(1-0.05)/12,
                     gamma2 = -log(1-0.05)/12,
                     accrualDuration = NA,
                     followupTime = 18, fixedFollowup = FALSE)
    testthat::expect_equal(
      round(l$resultsUnderH1$overallResults$accrualDuration, 2), 23.58)
  })


testthat::test_that(
  "follow-up time given power and accrual duration", {
    l = lrsamplesize(beta = 0.2, kMax = 2,
                     informationRates = c(0.8, 1),
                     alpha = 0.025, typeAlphaSpending = "sfOF",
                     accrualTime = seq(0, 8),
                     accrualIntensity = 26/9*seq(1, 9),
                     piecewiseSurvivalTime = c(0, 6),
                     stratumFraction = c(0.2, 0.8),
                     lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                     lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                     gamma1 = -log(1-0.05)/12,
                     gamma2 = -log(1-0.05)/12,
                     accrualDuration = 22,
                     followupTime = NA, fixedFollowup = FALSE)
    testthat::expect_equal(
      round(l$resultsUnderH1$overallResults$followupTime, 2), 21.55)
  })


testthat::test_that(
  "absolute accrual intensity given power, accrual duration, follow-up time,
  and relative accrual intensity", {
    l = lrsamplesize(beta = 0.2, kMax = 2,
                     informationRates = c(0.8, 1),
                     alpha = 0.025, typeAlphaSpending = "sfOF",
                     accrualTime = seq(0, 8),
                     accrualIntensity = 26/9*seq(1, 9),
                     piecewiseSurvivalTime = c(0, 6),
                     stratumFraction = c(0.2, 0.8),
                     lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
                     lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
                     gamma1 = -log(1-0.05)/12,
                     gamma2 = -log(1-0.05)/12,
                     accrualDuration = 22,
                     followupTime = 18, fixedFollowup = FALSE)
    testthat::expect_equal(
      round(l$resultsUnderH1$settings$accrualIntensity, 2),
      c(3.21, 6.42, 9.63, 12.84, 16.05, 19.26, 22.47, 25.68, 28.89))
  })

