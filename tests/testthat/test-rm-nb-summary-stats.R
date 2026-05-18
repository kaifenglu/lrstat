testthat::test_that("numeric results match hand-calculated examples", {
  rm1 <- rmstat(
    time = c(22, 40),
    milestone = 18,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  testthat::expect_equal(rm1$rmst1, c(10.8538733765985, 10.8538733765985), tolerance = 1e-6)
  testthat::expect_equal(rm1$rmst2, c(9.94809670616853, 9.94809670616853), tolerance = 1e-6)
  testthat::expect_equal(rm1$rmstDiff, c(0.905776670429968, 0.905776670429968), tolerance = 1e-6)
  testthat::expect_equal(rm1$rmstDiffZ, c(1.22254966257274, 1.466502154925), tolerance = 1e-6)

  nb1 <- nbstat(
    time = c(1, 1.25, 2, 3, 4),
    accrualIntensity = 1956 / 1.25,
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = 0.7 * 0.125,
    lambda2 = 0.125,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = 2.75
  )

  testthat::expect_equal(nb1$resultsUnderH1$rateRatio, rep(0.7, 5), tolerance = 1e-12)
  testthat::expect_equal(nb1$resultsUnderH1$zlogRR, c(-1.38715383992769, -1.68408165565445, -2.25139323856076, -2.60849129408142, -2.8062966871733), tolerance = 1e-6)
  testthat::expect_equal(nb1$resultsUnderH0$varianceRatio, rep(1, 5), tolerance = 1e-12)

  n1 <- 7
  n2 <- 8
  meanDiff <- 0.444
  stDev <- 1.201
  m <- n1 + n2 - 2
  ntilde <- n1 * n2 / (n1 + n2)
  tstat <- sqrt(ntilde) * meanDiff / stDev

  hg <- hedgesg(tstat, m, ntilde)
  testthat::expect_equal(hg$g, 0.347945299478477, tolerance = 1e-6)
  testthat::expect_equal(hg$varg, 0.272513478926849, tolerance = 1e-6)
  testthat::expect_equal(hg$lower, -0.675211266693981, tolerance = 1e-6)
  testthat::expect_equal(hg$upper, 1.37110186565094, tolerance = 1e-6)
})


testthat::test_that("rmst: Rd example regression and interval contracts", {
  # Rd example regression check.
  rm_rd <- rmst(
    t1 = 0,
    t2 = 7,
    piecewiseSurvivalTime = c(0, 6),
    lambda = c(0.0533, 0.0309)
  )

  expected_rd <- (1 - exp(-0.0533 * 6)) / 0.0533 +
    exp(-0.0533 * 6) * (1 - exp(-0.0309 * 1)) / 0.0309
  testthat::expect_equal(rm_rd, expected_rd, tolerance = 1e-10)

  # Additivity over adjacent intervals: RMST(0, 7) = RMST(0, 6) + RMST(6, 7).
  rm_0_6 <- rmst(t1 = 0, t2 = 6, piecewiseSurvivalTime = c(0, 6), lambda = c(0.0533, 0.0309))
  rm_6_7 <- rmst(t1 = 6, t2 = 7, piecewiseSurvivalTime = c(0, 6), lambda = c(0.0533, 0.0309))
  testthat::expect_equal(rm_rd, rm_0_6 + rm_6_7, tolerance = 1e-12)

  # Degenerate interval has zero restricted mean survival.
  rm_zero <- rmst(t1 = 5, t2 = 5, piecewiseSurvivalTime = c(0, 6), lambda = c(0.0533, 0.0309))
  testthat::expect_equal(rm_zero, 0, tolerance = 1e-12)
})


testthat::test_that("edge cases with small n and tied outcomes are handled", {
  # Tied outcome effect estimate gives t=0 and thus g=0.
  hg0 <- hedgesg(tstat = 0, m = 4, ntilde = 2)
  testthat::expect_equal(hg0$g, 0, tolerance = 1e-12)
  testthat::expect_true(hg0$lower < 0 && hg0$upper > 0)

  # Small df remains finite.
  hg_small <- hedgesg(tstat = 1.2, m = 2, ntilde = 1.5)
  testthat::expect_true(is.finite(hg_small$g))
  testthat::expect_true(is.finite(hg_small$varg))

  # Equal rates imply rate ratio 1 and near-zero z-statistics.
  nb_equal <- nbstat(
    time = c(1, 2),
    accrualIntensity = 200,
    kappa1 = 1,
    kappa2 = 1,
    lambda1 = 0.1,
    lambda2 = 0.1,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1,
    followupTime = 1
  )

  testthat::expect_equal(nb_equal$resultsUnderH1$rateRatio, c(1, 1), tolerance = 1e-12)
  testthat::expect_equal(nb_equal$resultsUnderH1$zlogRR, c(0, 0), tolerance = 1e-10)

  # rmstat contracts stay consistent for early/late cuts.
  rm_equal <- rmstat(
    time = c(7, 10),
    milestone = 6,
    allocationRatioPlanned = 1,
    accrualTime = 0,
    accrualIntensity = 30,
    piecewiseSurvivalTime = 0,
    stratumFraction = 1,
    lambda1 = 0.05,
    lambda2 = 0.05,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 5,
    followupTime = 6,
    fixedFollowup = FALSE
  )

  testthat::expect_equal(rm_equal$rmstDiff, c(0, 0), tolerance = 1e-6)
  testthat::expect_equal(rm_equal$rmstDiffZ, c(0, 0), tolerance = 1e-6)
})
