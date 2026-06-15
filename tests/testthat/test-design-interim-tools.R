testthat::test_that("getBound variants: documented numeric examples", {
  b <- getBound(
    k = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )
  bm <- getBound_multiarm(
    M = 2,
    k = 3,
    informationRates = seq(1, 3) / 3,
    alpha = 0.025,
    typeAlphaSpending = "OF"
  )
  bs <- getBound_seamless(
    M = 2,
    k = 2,
    informationRates = seq(1, 3) / 3,
    alpha = 0.025,
    typeAlphaSpending = "OF"
  )

  testthat::expect_equal(b, c(2.96258804272755, 1.96859564285759), tolerance = 1e-8)
  testthat::expect_equal(
    bm,
    c(3.88656166704350, 2.74821411026615, 2.24390742468965),
    tolerance = 1e-8
  )
  testthat::expect_equal(
    bs,
    c(3.77660487342257, 2.67046291585926, 2.18042384029337),
    tolerance = 1e-8
  )
})


testthat::test_that("getBound variants: reject invalid information rates", {
  testthat::expect_error(
    getBound(k = 2, informationRates = c(0.7, 0.6), alpha = 0.025),
    regexp = "increasing|strictly"
  )
  testthat::expect_error(
    getBound_multiarm(M = 2, k = 3, informationRates = c(0.3, 0.9, 1.2), alpha = 0.025),
    regexp = "<= 1|less than or equal|information"
  )
  testthat::expect_error(
    getBound_seamless(M = 2, k = 2, informationRates = c(0.6, 0.6, 1), alpha = 0.025),
    regexp = "increasing|strictly"
  )
})


testthat::test_that("exitprob: documented numeric example and structure", {
  out <- exitprob(
    b = c(3.471, 2.454, 2.004),
    theta = -log(0.6),
    I = c(50, 100, 150) / 4
  )

  testthat::expect_true(is.list(out))
  testthat::expect_named(out, c("exitProbUpper", "exitProbLower"))
  testthat::expect_equal(
    out$exitProbUpper,
    c(0.0479604991629151, 0.4927446216314498, 0.3327286203970084),
    tolerance = 1e-8
  )
  testthat::expect_equal(
    out$exitProbLower,
    c(5.30237470793893e-23, 2.42819712557121e-26, 0.126566482366414),
    tolerance = 1e-8
  )
  testthat::expect_true(all(cumsum(out$exitProbUpper) <= 1 + 1e-10))
  testthat::expect_true(all(cumsum(out$exitProbLower) <= 1 + 1e-10))
})


testthat::test_that("exitprob_multiarm: null and alternative examples", {
  I <- 95 / 2 * seq(1, 3) / 3
  b <- c(3.886562, 2.748214, 2.243907)

  p0 <- exitprob_multiarm(M = 2, theta = c(0, 0), kMax = 3, b = b, I = I)
  p1 <- exitprob_multiarm(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b, I = I)

  testthat::expect_named(p0, c("exitProbUpper", "exitProbLower"))
  testthat::expect_equal(
    p0$exitProbUpper,
    c(0.000100746545588226, 0.005707408864745522, 0.019191884452874675),
    tolerance = 1e-8
  )
  testthat::expect_equal(
    p1$exitProbUpper,
    c(0.0313048326907939, 0.5197405462908953, 0.3511726074563606),
    tolerance = 1e-8
  )
  testthat::expect_true(sum(p1$exitProbUpper) > sum(p0$exitProbUpper))
})


testthat::test_that("exitprob_seamless: documented examples and dimensions", {
  I <- 110 / 2 * seq(1, 3) / 3
  b <- c(3.776605, 2.670463, 2.180424)

  p0 <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = b, I = I)

  a <- c(0, 0.5, b[3])
  p1 <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2, b = b, a = a, I = I)

  testthat::expect_named(
    p0,
    c("exitProbUpper", "exitProbLower", "exitProbByArmUpper",
      "exitProbByArmLower", "selectAsBest")
  )
  testthat::expect_equal(dim(p0$exitProbByArmUpper), c(3, 2))
  testthat::expect_equal(dim(p0$exitProbByArmLower), c(3, 2))
  testthat::expect_length(p0$selectAsBest, 2)
  testthat::expect_equal(sum(p0$selectAsBest), 1, tolerance = 1e-8)

  testthat::expect_equal(
    p0$exitProbUpper,
    c(0.000157275641149868, 0.006485856535569953, 0.018356873822775486),
    tolerance = 1e-8
  )
  testthat::expect_equal(
    p1$exitProbLower,
    c(0.00780314438316704, 0.00623740048827642, 0.08795064352747015),
    tolerance = 1e-8
  )
})


testthat::test_that("getDurationFromNevents: structure and documented scenario", {
  out <- getDurationFromNevents(
    nevents = 80,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    fixedFollowup = FALSE
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_named(
    out,
    c("nevents", "fixedFollowup", "accrualDuration", "subjects", "followupTime", "studyDuration")
  )
  testthat::expect_equal(nrow(out), 23)
  testthat::expect_true(all(diff(out$accrualDuration) > 0))

  # Regression checks at first/middle/last grid points.
  testthat::expect_equal(out$accrualDuration[c(1, 12, 23)],
                         c(7.30797497709908, 11.69545054827838, 16.08292611945767),
                         tolerance = 1e-8)
  testthat::expect_equal(out$subjects[c(1, 12, 23)],
                         c(88.0065328040677, 200.0817142552378, 314.1560791058996),
                         tolerance = 1e-8)
})


testthat::test_that("getNeventsFromHazardRatio: documented scenario and monotonicity", {
  n1 <- getNeventsFromHazardRatio(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP",
    hazardRatio = 0.673
  )

  n_better <- getNeventsFromHazardRatio(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP",
    hazardRatio = 0.60
  )

  n_worse <- getNeventsFromHazardRatio(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP",
    hazardRatio = 0.80
  )

  testthat::expect_equal(n1, 234)
  testthat::expect_true(n_better < n1)
  testthat::expect_true(n_worse > n1)
})


testthat::test_that("getCI and getRCI: Rd example numeric regression", {
  n <- 120
  L <- 2
  n1 <- n * 2 / 3
  delta1 <- 7
  sigma1 <- 20
  zL <- delta1 / sqrt(4 / n1 * sigma1^2)

  ci <- getCI(
    L = L,
    zL = zL,
    IMax = n / (4 * sigma1^2),
    informationRates = c(1 / 3, 2 / 3),
    alpha = 0.05,
    typeAlphaSpending = "sfHSD",
    parameterAlphaSpending = -4
  )

  rci <- getRCI(
    L = L,
    zL = zL,
    IMax = n / (4 * sigma1^2),
    informationRates = c(1 / 3, 2 / 3),
    alpha = 0.05,
    typeAlphaSpending = "sfHSD",
    parameterAlphaSpending = -4
  )

  testthat::expect_s3_class(ci, "data.frame")
  testthat::expect_named(ci, c("pvalue", "thetahat", "cilevel", "lower", "upper"))
  testthat::expect_equal(round(ci$pvalue, 4), 0.0593)
  testthat::expect_equal(round(ci$thetahat, 4), 6.9882)
  testthat::expect_equal(ci$cilevel, 0.90)
  testthat::expect_equal(round(ci$lower, 4), -0.3764)
  testthat::expect_equal(round(ci$upper, 4), 14.3479)

  testthat::expect_s3_class(rci, "data.frame")
  testthat::expect_named(rci, c("pvalue", "thetahat", "cilevel", "lower", "upper"))
  testthat::expect_equal(round(rci$pvalue, 4), 0.2549)
  testthat::expect_equal(round(rci$thetahat, 4), 7.0000)
  testthat::expect_equal(rci$cilevel, 0.90)
  testthat::expect_equal(round(rci$lower, 4), -3.2367)
  testthat::expect_equal(round(rci$upper, 4), 17.2367)
  testthat::expect_true(ci$lower < ci$upper)
  testthat::expect_true(rci$lower < rci$upper)
})


testthat::test_that("getCI family: reject invalid information rates", {
  testthat::expect_error(
    getCI(L = 2, zL = 1.5, IMax = 10, informationRates = c(0.8, 0.7), alpha = 0.05),
    regexp = "increasing|strictly|information"
  )

  testthat::expect_error(
    getCI_multiarm(L = 2, zL = c(2.0, 2.1), M = 2, r = 1, corr_known = FALSE,
               IMax = 75, informationRates = c(0.7, 0.6), alpha = 0.025),
    regexp = "increasing|strictly|information"
  )

  testthat::expect_error(
    getCI_seamless(L = 2, zL = 2.0, M = 2, r = 1, corr_known = FALSE,
                   IMax = 75, informationRates = c(0.5, 0.5, 1), alpha = 0.025),
    regexp = "increasing|strictly|information"
  )
})


testthat::test_that("getCI_multiarm and getCI_seamless: Rd example numeric regression", {
  ci_multiarm <- getCI_multiarm(
    L = 2,
    zL = c(2.075, 2.264),
    M = 2,
    r = 1,
    corr_known = FALSE,
    IMax = 300 / 4,
    informationRates = c(1 / 2, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )

  ci_seam <- getCI_seamless(
    L = 2,
    zL = 2.075,
    M = 2,
    r = 1,
    corr_known = FALSE,
    IMax = 300 / 4,
    informationRates = c(1 / 3, 2 / 3, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(ci_multiarm, "data.frame")
  testthat::expect_named(ci_multiarm,
                         c("level", "index", "pvalue", "thetahat", "cilevel", "lower", "upper"))
  testthat::expect_equal(ci_multiarm$level, c(2L, 1L))
  testthat::expect_equal(ci_multiarm$index, c(2L, 1L))
  testthat::expect_equal(round(ci_multiarm$pvalue, 4), c(0.0241, 0.0196))
  testthat::expect_equal(round(ci_multiarm$thetahat, 4), c(0.1979, 0.2389))
  testthat::expect_equal(round(ci_multiarm$lower, 4), c(0.0017, 0.0119))
  testthat::expect_equal(round(ci_multiarm$upper, 4), c(0.4874, 0.4657))

  testthat::expect_s3_class(ci_seam, "data.frame")
  testthat::expect_named(ci_seam, c("pvalue", "thetahat", "cilevel", "lower", "upper"))
  testthat::expect_equal(round(ci_seam$pvalue, 4), 0.0342)
  testthat::expect_equal(round(ci_seam$thetahat, 4), 0.2006)
  testthat::expect_equal(round(ci_seam$lower, 4), -0.0153)
  testthat::expect_equal(round(ci_seam$upper, 4), 0.4652)
})


testthat::test_that("getCP family: Rd example numeric regression", {
  trialsdt <- as.Date("2020-03-04")
  iadt <- c(as.Date("2022-02-01"), as.Date("2022-11-01"))
  mo1 <- as.numeric(iadt - trialsdt + 1) / 30.4375

  N <- 521
  Ta <- 17.94
  Ta1 <- 8
  enrate <- N / (Ta - Ta1 / 2)

  lam1 <- log(2) / 16.7
  t1 <- 5
  hr <- 0.7
  lam2 <- hr * lam1
  gam <- -log(1 - 0.05) / 12

  mo2 <- caltime(
    nevents = c(298, 335),
    allocationRatioPlanned = 1,
    accrualTime = seq(0, Ta1),
    accrualIntensity = enrate * seq(1, Ta1 + 1) / (Ta1 + 1),
    piecewiseSurvivalTime = c(0, t1),
    lambda1 = c(lam1, lam2),
    lambda2 = c(lam1, lam1),
    gamma1 = gam,
    gamma2 = gam,
    accrualDuration = Ta,
    followupTime = 1000
  )

  lr1 <- lrstat(
    time = c(mo1, mo2),
    accrualTime = seq(0, Ta1),
    accrualIntensity = enrate * seq(1, Ta1 + 1) / (Ta1 + 1),
    piecewiseSurvivalTime = c(0, t1),
    lambda1 = c(lam1, lam2),
    lambda2 = c(lam1, lam1),
    gamma1 = gam,
    gamma2 = gam,
    accrualDuration = Ta,
    followupTime = 1000,
    predictTarget = 3
  )

  hr2 <- 0.81
  z2 <- (-log(hr2)) * sqrt(266 / 4)
  theta <- -log(lr1[, "HR"])

  cp <- getCP(
    INew = (335 - 266) / 4,
    L = 2,
    zL = z2,
    theta = theta,
    IMax = 298 / 4,
    kMax = 3,
    informationRates = c(179, 266, 298) / 298,
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )

  cp_multiarm <- getCP_multiarm(
    INew = 373 / 4,
    M = 2,
    r = 1,
    corr_known = FALSE,
    L = 1,
    zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
    theta = c(-log(0.91), -log(0.78)),
    IMax = 324 / 4,
    kMax = 2,
    informationRates = c(1 / 2, 1),
    alpha = 0.025,
    typeAlphaSpending = "OF",
    MNew = 1,
    selected = 2,
    rNew = 1
  )

  cp_seam <- getCP_seamless(
    INew = 198 / 4,
    M = 2,
    r = 1,
    corr_known = FALSE,
    L = 1,
    zL = -log(0.67) * sqrt(80 / 4),
    theta = -log(0.691),
    IMax = 120 / 4,
    K = 2,
    informationRates = c(1 / 3, 2 / 3, 1),
    alpha = 0.025,
    typeAlphaSpending = "OF",
    kNew = 1
  )

  testthat::expect_equal(round(cp, 4), c(0.3663, 0.5550))
  testthat::expect_equal(round(cp_multiarm, 4), c(0.4955, 0.8001))
  testthat::expect_equal(round(cp_seam, 4), c(0.4402, 0.6983))
  testthat::expect_true(cp[2] > cp[1])
  testthat::expect_true(cp_multiarm[2] > cp_multiarm[1])
  testthat::expect_true(cp_seam[2] > cp_seam[1])
})


testthat::test_that("getCP family: reject invalid information rates", {
  testthat::expect_error(
    getCP(INew = 1, L = 1, zL = 1.5, theta = c(0.1, 0.1), IMax = 10, kMax = 2,
          informationRates = c(1.1, 1), alpha = 0.025),
    regexp = "<= 1|increasing|information"
  )

  testthat::expect_error(
    getCP_multiarm(INew = 1, M = 2, r = 1, corr_known = FALSE, L = 1,
               zL = c(1.5, 1.6), theta = c(0.1, 0.2), IMax = 10, kMax = 2,
               informationRates = c(0.8, 0.7), alpha = 0.025),
    regexp = "increasing|strictly|information"
  )

  testthat::expect_error(
    getCP_seamless(INew = 1, M = 2, r = 1, corr_known = FALSE, L = 1,
                   zL = 1.5, theta = 0.1, IMax = 10, K = 2,
                   informationRates = c(0.6, 0.6, 1), alpha = 0.025),
    regexp = "increasing|strictly|information"
  )
})


testthat::test_that("getADCI, getADRCI, and adaptive variants: Rd example numeric regression", {
  delta <- 15
  sigma <- 50

  des1 <- getDesignMeanDiff(
    beta = 0.1,
    meanDiff = delta,
    stDev = sigma,
    kMax = 3,
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )

  s1 <- des1$byStageResults$informationRates
  n <- des1$overallResults$numberOfSubjects

  L <- 1
  nL <- des1$byStageResults$numberOfSubjects[L]
  deltahat <- 8
  sigmahat <- 55
  sedeltahat <- sigmahat * sqrt(4 / nL)
  zL <- deltahat / sedeltahat
  deltaNew <- 10

  des2 <- adaptDesign(
    betaNew = 0.1,
    L = L,
    zL = zL,
    theta = deltaNew,
    IMax = n / (4 * sigma^2),
    kMax = 3,
    informationRates = s1,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    MullerSchafer = TRUE,
    kNew = 2,
    typeAlphaSpendingNew = "sfP"
  )

  INew <- des2$secondaryTrial$maxInformation
  nNew <- ceiling(INew * 4 * sigma^2)
  s2 <- des2$secondaryTrial$informationRates

  Lc <- 2
  deltahatc <- 9.5
  sigmahatc <- 52.759
  L2 <- Lc - L
  nL2 <- nNew * s2[L2]
  nc <- nL + nL2
  sedeltahatc <- sigmahatc * sqrt(4 / nc)
  zLc <- deltahatc / sedeltahatc

  adc <- getADCI(
    L = L,
    zL = zL,
    IMax = n / (4 * sigmahatc^2),
    kMax = 3,
    informationRates = s1,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    MullerSchafer = TRUE,
    Lc = Lc,
    zLc = zLc,
    INew = nNew / (4 * sigmahatc^2),
    informationRatesNew = s2,
    typeAlphaSpendingNew = "sfP"
  )

  adrc <- getADRCI(
    L = L,
    zL = zL,
    IMax = n / (4 * sigmahatc^2),
    kMax = 3,
    informationRates = s1,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    MullerSchafer = TRUE,
    Lc = Lc,
    zLc = zLc,
    INew = nNew / (4 * sigmahatc^2),
    informationRatesNew = s2,
    typeAlphaSpendingNew = "sfP"
  )

  adc_multiarm <- getADCI_multiarm(
    M = 2,
    r = 1,
    corr_known = FALSE,
    L = 1,
    zL = c(2.075, 2.264),
    IMax = 300 / 4,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    MNew = 1,
    selected = 2,
    rNew = 1,
    Lc = 2,
    zLc = 1.667,
    INew = 374 / 4
  )

  adc_seam <- getADCI_seamless(
    M = 2,
    r = 1,
    corr_known = FALSE,
    L = 1,
    zL = -log(0.67) * sqrt(80 / 4),
    IMax = 120 / 4,
    K = 2,
    informationRates = c(1 / 3, 2 / 3, 1),
    alpha = 0.025,
    typeAlphaSpending = "OF",
    Lc = 2,
    zLc = -log(0.677) * sqrt(236 / 4),
    INew = 236 / 4
  )

  testthat::expect_s3_class(adc, "data.frame")
  testthat::expect_equal(round(adc$pvalue, 4), 0.0124)
  testthat::expect_equal(round(adc$thetahat, 4), 9.3496)
  testthat::expect_equal(round(adc$lower, 4), 1.2165)
  testthat::expect_equal(round(adc$upper, 4), 17.3653)

  testthat::expect_s3_class(adrc, "data.frame")
  testthat::expect_equal(round(adrc$pvalue, 4), 0.0167)
  testthat::expect_equal(round(adrc$thetahat, 4), 10.0537)
  testthat::expect_equal(round(adrc$lower, 4), 0.6303)
  testthat::expect_equal(round(adrc$upper, 4), 18.1808)

  testthat::expect_s3_class(adc_multiarm, "data.frame")
  testthat::expect_equal(adc_multiarm$level, 1L)
  testthat::expect_equal(adc_multiarm$index, 2L)
  testthat::expect_equal(round(adc_multiarm$pvalue, 4), 0.0255)
  testthat::expect_equal(round(adc_multiarm$thetahat, 4), 0.1713)
  testthat::expect_equal(round(adc_multiarm$lower, 4), -0.0008)
  testthat::expect_equal(round(adc_multiarm$upper, 4), 0.3530)

  testthat::expect_s3_class(adc_seam, "data.frame")
  testthat::expect_equal(round(adc_seam$pvalue, 4), 0.0085)
  testthat::expect_equal(round(adc_seam$thetahat, 4), 0.2557)
  testthat::expect_equal(round(adc_seam$lower, 4), 0.0557)
  testthat::expect_equal(round(adc_seam$upper, 4), 0.4421)
})


testthat::test_that("getADCI family: reject malformed new information rates", {
  testthat::expect_error(
    getADCI(L = 1, zL = 1.5, IMax = 10, kMax = 2, informationRates = c(0.8, 0.7),
            alpha = 0.025, Lc = 2, zLc = 1.7, INew = 3),
    regexp = "increasing|strictly|information"
  )

  testthat::expect_error(
    getADCI_multiarm(M = 2, r = 1, corr_known = FALSE, L = 1, zL = c(1.5, 1.6),
                 IMax = 10, kMax = 2, informationRates = c(0.5, 1), alpha = 0.025,
                 MNew = 1, selected = 3, rNew = 1, Lc = 2, zLc = 1.7, INew = 3),
    regexp = "Invalid value in selected|selected"
  )

  testthat::expect_error(
    getADCI_seamless(M = 2, r = 1, corr_known = FALSE, L = 1, zL = 1.5,
                     IMax = 10, K = 2, informationRates = c(0.5, 0.4, 1),
                     alpha = 0.025, Lc = 2, zLc = 1.7, INew = 3),
    regexp = "increasing|strictly|information"
  )
})
