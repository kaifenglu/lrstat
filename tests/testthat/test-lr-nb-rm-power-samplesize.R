testthat::test_that("LR equivalence and Schoenfeld outputs match documented examples", {
  lr_eq <- lrpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    hazardRatioLower = 0.71,
    hazardRatioUpper = 1.4,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 100 / 9 * seq(1, 9),
    lambda1 = 0.0533,
    lambda2 = 0.0533,
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  lr_ss_eq <- lrsamplesizeequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    hazardRatioLower = 0.71,
    hazardRatioUpper = 1.4,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0533),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = NA,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  lr_sch <- lrschoenfeld(
    beta = 0.1,
    kMax = 2,
    alpha = 0.025,
    hazardRatioH0 = 1,
    allocationRatioPlanned = 1,
    accrualIntensity = 20,
    hazardRatio = 0.3,
    lambda2 = 1.9 / 12,
    gamma1 = -log(1 - 0.1) / 24,
    gamma2 = -log(1 - 0.1) / 24,
    fixedFollowup = 0,
    rounding = 1,
    calibrate = 0,
    maxNumberOfIterations = 1000,
    seed = 12345
  )

  testthat::expect_s3_class(lr_eq, "lrpowerequiv")
  testthat::expect_equal(round(lr_eq$overallResults$overallReject, 6), 0.999992)
  testthat::expect_equal(lr_eq$overallResults$numberOfSubjects, 1800)

  testthat::expect_s3_class(lr_ss_eq, "lrpowerequiv")
  testthat::expect_equal(round(lr_ss_eq$overallResults$overallReject, 6), 0.800927)
  testthat::expect_equal(lr_ss_eq$overallResults$numberOfSubjects, 420)
  testthat::expect_equal(round(lr_ss_eq$overallResults$accrualDuration, 5), 20.15385)

  testthat::expect_named(lr_sch, c("analyticalResults", "simulationResults"))
  testthat::expect_equal(round(lr_sch$analyticalResults$overallResults$overallReject, 6), 0.908518)
  testthat::expect_equal(lr_sch$analyticalResults$overallResults$numberOfEvents, 30)
  testthat::expect_equal(lr_sch$analyticalResults$overallResults$numberOfSubjects, 78)
})


testthat::test_that("NB power and sample-size outputs match documented examples", {
  nb <- nbpower(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    accrualIntensity = 1956 / 1.25,
    stratumFraction = c(0.2, 0.8),
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = c(0.7 * 0.125, 0.75 * 0.25),
    lambda2 = c(0.125, 0.25),
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = 2.75,
    fixedFollowup = FALSE,
    nullVariance = 1
  )

  nb1 <- nbpower1s(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    lambdaH0 = 0.125,
    accrualIntensity = 500,
    stratumFraction = c(0.2, 0.8),
    kappa = c(3, 5),
    lambda = c(0.0875, 0.085),
    gamma = 0,
    accrualDuration = 1.25,
    followupTime = 2.75,
    fixedFollowup = FALSE
  )

  nb_eq <- nbpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    rateRatioLower = 2 / 3,
    rateRatioUpper = 3 / 2,
    accrualIntensity = 1956 / 1.25,
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = 0.125,
    lambda2 = 0.125,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = 2.75,
    fixedFollowup = FALSE
  )

  nb_ss <- nbsamplesize(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    accrualIntensity = 1956 / 1.25,
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = 0.0875,
    lambda2 = 0.125,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = NA,
    fixedFollowup = FALSE
  )

  nb_ss1 <- nbsamplesize1s(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    lambdaH0 = 0.125,
    accrualIntensity = 500,
    stratumFraction = c(0.2, 0.8),
    kappa = c(3, 5),
    lambda = c(0.0875, 0.085),
    gamma = 0,
    accrualDuration = 1.25,
    followupTime = NA,
    fixedFollowup = FALSE
  )

  nb_sseq <- nbsamplesizeequiv(
    beta = 0.1,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    rateRatioLower = 2 / 3,
    rateRatioUpper = 3 / 2,
    accrualIntensity = 1956 / 1.25,
    stratumFraction = c(0.2, 0.8),
    kappa1 = c(3, 5),
    kappa2 = c(2, 3),
    lambda1 = c(0.125, 0.165),
    lambda2 = c(0.135, 0.175),
    gamma1 = -log(1 - 0.05),
    gamma2 = -log(1 - 0.10),
    accrualDuration = 1.25,
    followupTime = NA,
    fixedFollowup = FALSE
  )

  testthat::expect_equal(round(nb$overallResults$overallReject, 6), 0.731577)
  testthat::expect_equal(nb$overallResults$numberOfSubjects, 1956)

  testthat::expect_equal(round(nb1$overallResults$overallReject, 6), 0.915213)
  testthat::expect_equal(nb1$overallResults$numberOfSubjects, 625)

  testthat::expect_equal(round(nb_eq$overallResults$overallReject, 6), 0.89954)
  testthat::expect_equal(nb_eq$overallResults$numberOfSubjects, 1956)

  testthat::expect_equal(round(nb_ss$resultsUnderH1$overallResults$overallReject, 4), 0.8)
  testthat::expect_equal(round(nb_ss$resultsUnderH1$overallResults$followupTime, 6), 2.753248)
  testthat::expect_equal(nb_ss$resultsUnderH1$overallResults$numberOfSubjects, 1956)

  testthat::expect_equal(round(nb_ss1$resultsUnderH1$overallResults$overallReject, 4), 0.8)
  testthat::expect_equal(round(nb_ss1$resultsUnderH1$overallResults$followupTime, 6), 1.113623)
  testthat::expect_equal(nb_ss1$resultsUnderH1$overallResults$numberOfSubjects, 625)

  testthat::expect_equal(round(nb_sseq$overallResults$overallReject, 4), 0.9)
  testthat::expect_equal(round(nb_sseq$overallResults$followupTime, 6), 2.082765)
  testthat::expect_equal(nb_sseq$overallResults$numberOfSubjects, 1956)
})


testthat::test_that("RM power and sample-size outputs match documented examples", {
  rm <- rmpower(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
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

  rm1 <- rmpower1s(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstH0 = 10,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    gamma = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  rm_eq <- rmpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstDiffLower = -2,
    rmstDiffUpper = 2,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 29 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  rm_ss <- rmsamplesize(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 100 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = NA,
    fixedFollowup = FALSE
  )

  rm_ss1 <- rmsamplesize1s(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstH0 = 10,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    gamma = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = NA,
    fixedFollowup = FALSE
  )

  rm_sseq <- rmsamplesizeequiv(
    beta = 0.1,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstDiffLower = -2,
    rmstDiffUpper = 2,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = NA,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  testthat::expect_equal(round(rm$overallResults$overallReject, 6), 0.30423)
  testthat::expect_equal(rm$overallResults$numberOfSubjects, 468)

  testthat::expect_equal(round(rm1$overallResults$overallReject, 6), 0.755621)
  testthat::expect_equal(rm1$overallResults$numberOfSubjects, 468)

  testthat::expect_equal(round(rm_eq$overallResults$overallReject, 6), 0.939668)
  testthat::expect_equal(rm_eq$overallResults$numberOfSubjects, 522)
  testthat::expect_equal(round(rm_eq$overallResults$rmstDiff, 6), 0)

  testthat::expect_equal(round(rm_ss$resultsUnderH1$overallResults$overallReject, 4), 0.8)
  testthat::expect_equal(round(rm_ss$resultsUnderH1$overallResults$followupTime, 6), 7.488332)
  testthat::expect_equal(rm_ss$resultsUnderH1$overallResults$numberOfSubjects, 1800)

  testthat::expect_equal(round(rm_ss1$resultsUnderH1$overallResults$overallReject, 4), 0.8)
  testthat::expect_equal(round(rm_ss1$resultsUnderH1$overallResults$followupTime, 5), 15.15319)
  testthat::expect_equal(rm_ss1$resultsUnderH1$overallResults$numberOfSubjects, 522)

  testthat::expect_equal(round(rm_sseq$overallResults$overallReject, 4), 0.9)
  testthat::expect_equal(round(rm_sseq$overallResults$accrualDuration, 5), 21.53846)
  testthat::expect_equal(rm_sseq$overallResults$numberOfSubjects, 456)
})


testthat::test_that("equivalence margins behave correctly at null boundaries", {
  lr_inside <- lrpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    hazardRatioLower = 0.71,
    hazardRatioUpper = 1.4,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 100 / 9 * seq(1, 9),
    lambda1 = 0.0533,
    lambda2 = 0.0533,
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  lr_boundary <- lrpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    hazardRatioLower = 0.71,
    hazardRatioUpper = 1.4,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 100 / 9 * seq(1, 9),
    lambda1 = 1.4 * 0.0533,
    lambda2 = 0.0533,
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  nb_inside <- nbpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    rateRatioLower = 2 / 3,
    rateRatioUpper = 3 / 2,
    accrualIntensity = 1956 / 1.25,
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = 0.125,
    lambda2 = 0.125,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = 2.75,
    fixedFollowup = FALSE
  )

  nb_boundary <- nbpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    rateRatioLower = 2 / 3,
    rateRatioUpper = 3 / 2,
    accrualIntensity = 1956 / 1.25,
    kappa1 = 5,
    kappa2 = 5,
    lambda1 = 1.5 * 0.125,
    lambda2 = 0.125,
    gamma1 = 0,
    gamma2 = 0,
    accrualDuration = 1.25,
    followupTime = 2.75,
    fixedFollowup = FALSE
  )

  rm_inside <- rmpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstDiffLower = -2,
    rmstDiffUpper = 2,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 29 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  rm_boundary <- rmpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    rmstDiffLower = -2,
    rmstDiffUpper = 2,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 29 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    lambda2 = c(0.0533 / 1.4, 0.0533 / 1.4, 1.5 * 0.0533 / 1.4, 1.5 * 0.0533 / 1.4),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  testthat::expect_gt(lr_inside$overallResults$overallReject, lr_boundary$overallResults$overallReject)
  testthat::expect_equal(round(lr_boundary$overallResults$overallReject, 2), 0.05)

  testthat::expect_gt(nb_inside$overallResults$overallReject, nb_boundary$overallResults$overallReject)
  testthat::expect_equal(round(nb_boundary$overallResults$overallReject, 2), 0.05)

  testthat::expect_gt(rm_inside$overallResults$overallReject, rm_boundary$overallResults$overallReject)
  testthat::expect_true(rm_boundary$overallResults$overallReject < 0.2)
})


testthat::test_that("invalid hazard rates or follow-up values raise errors", {
  testthat::expect_error(
    lrpowerequiv(
      kMax = 2,
      informationRates = c(0.5, 1),
      alpha = 0.05,
      typeAlphaSpending = "sfOF",
      hazardRatioLower = 0.71,
      hazardRatioUpper = 1.4,
      allocationRatioPlanned = 1,
      accrualTime = seq(0, 8),
      accrualIntensity = 100 / 9 * seq(1, 9),
      lambda1 = -0.0533,
      lambda2 = 0.0533,
      gamma1 = -log(1 - 0.05) / 12,
      gamma2 = -log(1 - 0.05) / 12,
      accrualDuration = 22,
      followupTime = 18,
      fixedFollowup = FALSE
    )
  )

  testthat::expect_error(
    nbpower(
      kMax = 2,
      informationRates = c(0.5, 1),
      alpha = 0.025,
      typeAlphaSpending = "sfOF",
      accrualIntensity = 1956 / 1.25,
      stratumFraction = c(0.2, 0.8),
      kappa1 = 5,
      kappa2 = 5,
      lambda1 = c(-0.7 * 0.125, 0.75 * 0.25),
      lambda2 = c(0.125, 0.25),
      gamma1 = 0,
      gamma2 = 0,
      accrualDuration = 1.25,
      followupTime = 2.75,
      fixedFollowup = FALSE,
      nullVariance = 1
    )
  )

  testthat::expect_error(
    rmpower(
      kMax = 2,
      informationRates = c(0.8, 1),
      alpha = 0.025,
      typeAlphaSpending = "sfOF",
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
      followupTime = -1,
      fixedFollowup = FALSE
    )
  )
})



testthat::test_that("lrpower: power with alpha-spending, weighted", {
  l = lrpower(kMax = 2, informationRates = c(0.8, 1),
              alpha = 0.025, typeAlphaSpending = "sfOF",
              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
              accrualIntensity = 26/9*seq(1, 9),
              piecewiseSurvivalTime = c(0, 6),
              stratumFraction = c(0.2, 0.8),
              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
              gamma1 = -log(1-0.05)/12,
              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
              followupTime = 18, fixedFollowup = FALSE,
              rho1 = 0, rho2 = 1)
  testthat::expect_equal(round(l$overallResults$overallReject, 4), 0.9313)
})


testthat::test_that("lrpower: power for stratified analysis", {
  p1 = c(0.28, 0.13, 0.25, 0.34)
  p2 = c(0.28, 0.72)
  p3 = c(0.43, 0.37, 0.2)
  stratumFraction = p1 %x% p2 %x% p3
  stratumFraction = stratumFraction/sum(stratumFraction)
  theta1 = c(1, 2.127, 0.528, 0.413)
  theta2 = c(1, 0.438)
  theta3 = c(1, 0.614, 0.159)
  lambda2 = 0.009211*exp(log(theta1) %x% log(theta2) %x% log(theta3))
  caltime(nevents = 66, accrualDuration = 24, accrualIntensity = 12,
          stratumFraction = stratumFraction,
          lambda1 = 0.4466*lambda2, lambda2 = lambda2,
          followupTime = 100)
  l = lrpower(kMax = 3,
              informationRates = (1:3)/3,
              alpha = 0.025, typeAlphaSpending = "sfOF",
              accrualIntensity = 12,
              stratumFraction = stratumFraction,
              lambda1 = 0.4466*lambda2,
              lambda2 = lambda2,
              accrualDuration = 24,
              followupTime = 30.92,
              typeOfComputation = "schoenfeld")
  testthat::expect_equal(round(l$overallResults$overallReject, 4), 0.9024)
})



testthat::test_that("lrsamplesize: accrual duration given power and follow-up time", {
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


testthat::test_that("lrsamplesize: follow-up time given power and accrual duration", {
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


testthat::test_that("lrsamplesize: absolute accrual intensity given power, accrual duration, follow-up time, and relative accrual intensity", {
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



testthat::test_that("rmest: valid milestone time works", {
  library(dplyr, warn.conflicts = FALSE)
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time",
               event="status", milestone=161)

  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  rmst = 52.64545
  stderr = 19.8286

  testthat::expect_equal(round(df1$rmst, 5), rmst)
  testthat::expect_equal(round(df1$stderr, 4), stderr)
})


testthat::test_that("rmest: bias correction for variance estimate", {
  library(dplyr, warn.conflicts = FALSE)
  df1 <- rmest(aml %>% filter(x == "Maintained"), time="time",
               event="status", milestone=161, biascorrection=TRUE)

  # the following estimates can be obtained from SAS PROC LIFETEST
  #   proc lifetest data=aml rmst(bc tau=161);
  #     where x = "Maintained";
  #     time time*status(0);
  #   run;
  stderr = 21.4173

  testthat::expect_equal(round(df1$stderr, 4), stderr)
})
