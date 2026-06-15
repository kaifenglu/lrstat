testthat::test_that("getDesign: documented example returns stable design object", {
  design1 <- getDesign(
    beta = 0.2,
    theta = -log(0.7),
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP"
  )

  testthat::expect_s3_class(design1, "design")
  testthat::expect_named(design1, c("byStageResults", "overallResults", "settings"))

  # Printed summary reports power and max information; verify numerically.
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8000)
  testthat::expect_equal(round(design1$overallResults$information, 2), 71.97)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 2.963)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[2], 3), 1.969)
})


testthat::test_that("getDesign_multiarm: documented example returns stable multiarm object", {
  design1 <- getDesign_multiarm(
    beta = 0.1,
    theta = c(0.3, 0.5),
    M = 2,
    r = 1.0,
    kMax = 3,
    informationRates = seq(1, 3) / 3,
    alpha = 0.025,
    typeAlphaSpending = "OF"
  )

  testthat::expect_s3_class(design1, "multiarm")
  testthat::expect_named(design1,
                         c("byStageResults", "byLevelBounds", "overallResults", "settings"))

  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9000)
  testthat::expect_equal(round(design1$overallResults$information, 2), 47.15)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.887, 2.748, 2.244)
  )
})


testthat::test_that("getDesign_seamless: object contract and boundary dimensions", {
  design1 <- getDesign_seamless(
    beta = 0.1,
    theta = c(0.3, 0.5),
    M = 2,
    r = 1.0,
    K = 2,
    informationRates = seq(1, 3) / 3,
    alpha = 0.025,
    typeAlphaSpending = "OF"
  )

  testthat::expect_s3_class(design1, "seamless")
  testthat::expect_true("overallResults" %in% names(design1))
  testthat::expect_true("byStageResults" %in% names(design1))
  testthat::expect_true("settings" %in% names(design1))
  testthat::expect_equal(nrow(design1$byStageResults), 3)
  testthat::expect_true(all(design1$byStageResults$efficacyBounds > 0))
})


testthat::test_that("getDesign wrappers: invalid alpha is rejected", {
  testthat::expect_error(
    getDesign(beta = 0.2, theta = 0.4, kMax = 2, informationRates = c(0.5, 1), alpha = -0.1),
    regexp = "alpha|significance"
  )

  testthat::expect_error(
    getDesign_multiarm(beta = 0.1, theta = c(0.3, 0.5), M = 2, kMax = 3,
                   informationRates = seq(1, 3) / 3, alpha = -0.1),
    regexp = "alpha|significance"
  )

  testthat::expect_error(
    getDesign_seamless(beta = 0.1, theta = c(0.3, 0.5), M = 2, K = 2,
                       informationRates = seq(1, 3) / 3, alpha = -0.1),
    regexp = "alpha|significance"
  )
})


testthat::test_that("getDesignRiskDiff: Rd example numeric regression", {
  design1 <- getDesignRiskDiff(
    beta = 0.2,
    n = NA,
    pi1 = 0.1,
    pi2 = 0.15,
    kMax = 3,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    nullVariance = FALSE
  )

  testthat::expect_s3_class(design1, "designRiskDiff")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8002)
  testthat::expect_equal(round(design1$overallResults$numberOfSubjects, 0), 1384)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.712, 2.511, 1.993)
  )
})


testthat::test_that("getDesignOneProportion: Rd examples numeric regression", {
  design1 <- getDesignOneProportion(
    beta = 0.2,
    n = NA,
    piH0 = 0.15,
    pi = 0.25,
    kMax = 3,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignOneProportion(
    beta = 0.2,
    n = NA,
    piH0 = 0.15,
    pi = 0.25,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designOneProportion")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8012)
  testthat::expect_equal(round(design1$overallResults$numberOfSubjects, 0), 91)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[3], 3), 1.696)

  testthat::expect_s3_class(design2, "designOneProportion")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.8097)
  testthat::expect_equal(round(design2$overallResults$numberOfSubjects, 0), 110)
})


testthat::test_that("getDesignFisherExact: Rd example numeric regression", {
  design1 <- getDesignFisherExact(
    beta = 0.2,
    pi1 = 0.5,
    pi2 = 0.2,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(round(design1$power, 4), 0.8168)
  testthat::expect_equal(design1$n, 87)
})


testthat::test_that("getDesignMeanDiff: Rd examples numeric regression", {
  design1 <- getDesignMeanDiff(
    beta = NA,
    n = 456,
    meanDiff = 9,
    stDev = 32,
    kMax = 5,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP"
  )

  design2 <- getDesignMeanDiff(
    beta = 0.1,
    n = NA,
    meanDiff = 0.3,
    stDev = 1,
    normalApproximation = FALSE,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designMeanDiff")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.7421)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 456)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 4.883)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[5], 3), 2.031)

  testthat::expect_s3_class(design2, "designMeanDiff")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9000)
  testthat::expect_equal(round(design2$overallResults$numberOfSubjects, 0), 469)
})


testthat::test_that("getDesignWilcoxon: Rd examples numeric regression", {
  p_larger <- pnorm((8 - 2) / sqrt(2 * 25^2))

  design1 <- getDesignWilcoxon(
    beta = 0.1,
    n = NA,
    pLarger = p_larger,
    alpha = 0.025
  )

  design2 <- getDesignWilcoxon(
    beta = 0.1,
    n = NA,
    pLarger = p_larger,
    alpha = 0.025,
    kMax = 3,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designWilcoxon")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9002)
  testthat::expect_equal(round(design1$overallResults$numberOfSubjects, 0), 772)

  testthat::expect_s3_class(design2, "designWilcoxon")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9001)
  testthat::expect_equal(round(design2$overallResults$numberOfSubjects, 0), 781)
  testthat::expect_equal(
    round(design2$byStageResults$efficacyBounds, 3),
    c(3.713, 2.510, 1.993)
  )
})


testthat::test_that("getDesignOddsRatio: Rd example numeric regression", {
  design1 <- getDesignOddsRatio(
    beta = 0.1,
    n = NA,
    pi1 = 0.5,
    pi2 = 0.3,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designOddsRatio")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9012)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 210)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 1.645)
})


testthat::test_that("getDesignRiskRatio: Rd example numeric regression", {
  design1 <- getDesignRiskRatio(
    beta = 0.1,
    n = NA,
    pi1 = 0.5,
    pi2 = 0.3,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designRiskRatio")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9008)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 207)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 1.645)
})


testthat::test_that("getDesignOneMean: Rd examples numeric regression", {
  design1 <- getDesignOneMean(
    beta = 0.1,
    n = NA,
    meanH0 = 7,
    mean = 6,
    stDev = 2.5,
    kMax = 5,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    typeBetaSpending = "sfP"
  )

  design2 <- getDesignOneMean(
    beta = 0.1,
    n = NA,
    meanH0 = 7,
    mean = 6,
    stDev = 2.5,
    normalApproximation = FALSE,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designOneMean")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9016)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 85)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 4.877)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[5], 3), 2.031)

  testthat::expect_s3_class(design2, "designOneMean")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9016)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 68)
  testthat::expect_equal(round(design2$byStageResults$efficacyBounds[1], 3), 1.996)
})


testthat::test_that("getDesignANOVA: Rd example numeric regression", {
  design1 <- getDesignANOVA(
    beta = 0.1,
    ngroups = 4,
    means = c(1.5, 2.5, 2, 0),
    stDev = 3.5,
    allocationRatioPlanned = c(2, 2, 2, 1),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designANOVA")
  testthat::expect_equal(round(design1$power, 4), 0.9008)
  testthat::expect_equal(design1$n, 279)
  testthat::expect_equal(round(design1$effectsize, 4), 0.0516)
})


testthat::test_that("getDesignTwoOrdinal: Rd example numeric regression", {
  design1 <- getDesignTwoOrdinal(
    beta = 0.1,
    ncats = 4,
    pi1 = c(0.55, 0.3, 0.1),
    pi2 = c(0.214, 0.344, 0.251),
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designTwoOrdinal")
  testthat::expect_equal(round(design1$power, 4), 0.9030)
  testthat::expect_equal(design1$n, 67)
  testthat::expect_equal(round(design1$meanscore1, 4), 26.4055)
  testthat::expect_equal(round(design1$meanscore2, 4), 40.5945)
})


testthat::test_that("getDesignOddsRatioEquiv: Rd example numeric regression", {
  design1 <- getDesignOddsRatioEquiv(
    beta = 0.2,
    n = NA,
    oddsRatioLower = 0.8,
    oddsRatioUpper = 1.25,
    pi1 = 0.12,
    pi2 = 0.12,
    kMax = 3,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designOddsRatioEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8000)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 6637)
  testthat::expect_equal(round(design1$overallResults$information, 4), 175.2168)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.200, 2.141, 1.695)
  )
})


testthat::test_that("getDesignRiskRatioEquiv: Rd example numeric regression", {
  design1 <- getDesignRiskRatioEquiv(
    beta = 0.2,
    n = NA,
    riskRatioLower = 0.8,
    riskRatioUpper = 1.25,
    pi1 = 0.12,
    pi2 = 0.12,
    kMax = 3,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designRiskRatioEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8000)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 5140)
  testthat::expect_equal(round(design1$overallResults$information, 4), 175.2273)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.200, 2.141, 1.695)
  )
})


testthat::test_that("getDesignOneRateExact: Rd examples numeric regression", {
  design1 <- getDesignOneRateExact(
    n = 525,
    lambdaH0 = 0.049,
    lambda = 0.012,
    D = 0.5,
    alpha = 0.025
  )

  design2 <- getDesignOneRateExact(
    beta = 0.2,
    lambdaH0 = 0.2,
    lambda = 0.3,
    D = 1,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(round(design1$power, 4), 0.9002)
  testthat::expect_equal(round(design1$attainedAlpha, 4), 0.0117)
  testthat::expect_equal(design1$n, 525)
  testthat::expect_equal(design1$r, 5)

  testthat::expect_s3_class(design2, "data.frame")
  testthat::expect_equal(round(design2$power, 4), 0.8078)
  testthat::expect_equal(round(design2$attainedAlpha, 4), 0.0427)
  testthat::expect_equal(design2$n, 162)
  testthat::expect_equal(design2$r, 43)
})


testthat::test_that("getDesignANOVAContrast: Rd example numeric regression", {
  design1 <- getDesignANOVAContrast(
    beta = 0.1,
    ngroups = 4,
    means = c(1.5, 2.5, 2, 0),
    stDev = 3.5,
    contrast = c(1, 1, 1, -3),
    allocationRatioPlanned = c(2, 2, 2, 1),
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designANOVAContrast")
  testthat::expect_equal(round(design1$power, 4), 0.9003)
  testthat::expect_equal(design1$n, 265)
  testthat::expect_equal(round(design1$effectsize, 4), 0.0400)
  testthat::expect_equal(design1$meanContrast, 6)
})


testthat::test_that("getDesignRepeatedANOVA: Rd example numeric regression", {
  design1 <- getDesignRepeatedANOVA(
    beta = 0.1,
    ngroups = 4,
    means = c(1.5, 2.5, 2, 0),
    stDev = 5,
    corr = 0.2,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designRepeatedANOVA")
  testthat::expect_equal(round(design1$power, 4), 0.9027)
  testthat::expect_equal(design1$n, 83)
  testthat::expect_equal(round(design1$effectsize, 3), 0.175)
})


testthat::test_that("getDesignMeanDiffEquiv: Rd examples numeric regression", {
  design1 <- getDesignMeanDiffEquiv(
    beta = 0.1,
    n = NA,
    meanDiffLower = -1.3,
    meanDiffUpper = 1.3,
    meanDiff = 0,
    stDev = 2.2,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignMeanDiffEquiv(
    beta = 0.1,
    n = NA,
    meanDiffLower = -1.3,
    meanDiffUpper = 1.3,
    meanDiff = 0,
    stDev = 2.2,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designMeanDiffEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9024)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 128)
  testthat::expect_equal(round(design1$overallResults$information, 4), 6.6116)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.750, 2.540, 2.016, 1.720)
  )

  testthat::expect_s3_class(design2, "designMeanDiffEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9018)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 126)
  testthat::expect_equal(round(design2$byStageResults$efficacyBounds[1], 3), 1.657)
})


testthat::test_that("getDesignAgreement: Rd example numeric regression", {
  design1 <- getDesignAgreement(
    beta = 0.2,
    n = NA,
    ncats = 4,
    kappaH0 = 0.4,
    kappa = 0.6,
    p1 = c(0.1, 0.2, 0.3, 0.4),
    p2 = c(0.15, 0.2, 0.24, 0.41),
    rounding = TRUE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designAgreement")
  testthat::expect_equal(round(design1$power, 4), 0.8006)
  testthat::expect_equal(design1$n, 82)
})


testthat::test_that("getDesignEquiv: Rd example numeric regression", {
  design1 <- getDesignEquiv(
    beta = 0.2,
    thetaLower = log(0.8),
    thetaUpper = log(1.25),
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8000)
  testthat::expect_equal(round(design1$overallResults$information, 4), 173.2257)
})


testthat::test_that("getDesignLogistic: Rd example numeric regression", {
  x1 <- c(5, 10, 15, 20)
  px1 <- c(0.2, 0.3, 0.3, 0.2)
  x2 <- c(2, 4, 6)
  px2 <- c(0.4, 0.4, 0.2)
  nbins <- 10
  x3 <- qnorm(((1:nbins) - 0.5) / nbins) * 2 + 4
  px3 <- rep(1 / nbins, nbins)

  nconfigs <- length(x1) * length(x2) * length(x3)
  x <- expand.grid(x3 = x3, x2 = x2, x1 = x1)
  x <- as.matrix(x[, ncol(x):1])
  pconfigs <- as.numeric(px1 %x% px2 %x% px3)

  design1 <- getDesignLogistic(
    beta = 0.1,
    ncovariates = 3,
    nconfigs = nconfigs,
    x = x,
    pconfigs = pconfigs,
    oddsratios = c(1.2^(1/5), 1.4, 1.3),
    responseprob = 0.25,
    alpha = 0.1
  )

  testthat::expect_s3_class(design1, "designLogistic")
  testthat::expect_equal(round(design1$power, 4), 0.9002)
  testthat::expect_equal(design1$n, 1369)
  testthat::expect_equal(round(design1$effectsize, 4), 0.0063)
})


testthat::test_that("getDesignMeanRatio: Rd example numeric regression", {
  design1 <- getDesignMeanRatio(
    beta = 0.1,
    n = NA,
    meanRatio = 1.25,
    CV = 0.35,
    kMax = 3,
    alpha = 0.025,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designMeanRatio")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9009)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 99)
})


testthat::test_that("getDesignOneMultinom: Rd example numeric regression", {
  design1 <- getDesignOneMultinom(
    beta = 0.1,
    ncats = 3,
    piH0 = c(0.25, 0.25),
    pi = c(0.3, 0.4),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designOneMultinom")
  testthat::expect_equal(round(design1$power, 4), 0.9030)
  testthat::expect_equal(design1$n, 71)
})


testthat::test_that("getDesignMeanRatioEquiv: Rd examples numeric regression", {
  design1 <- getDesignMeanRatioEquiv(
    beta = 0.1,
    n = NA,
    meanRatioLower = 0.8,
    meanRatioUpper = 1.25,
    meanRatio = 1,
    CV = 0.35,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignMeanRatioEquiv(
    beta = 0.1,
    n = NA,
    meanRatioLower = 0.8,
    meanRatioUpper = 1.25,
    meanRatio = 1,
    CV = 0.35,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designMeanRatioEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9033)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 104)
  testthat::expect_equal(round(design1$overallResults$information, 4), 224.9946)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.750, 2.540, 2.016, 1.720)
  )

  testthat::expect_s3_class(design2, "designMeanRatioEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9005)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 102)
  testthat::expect_equal(round(design2$byStageResults$efficacyBounds[1], 3), 1.660)
})


testthat::test_that("getDesignRiskDiffEquiv: Rd example numeric regression", {
  design1 <- getDesignRiskDiffEquiv(
    beta = 0.2,
    n = NA,
    riskDiffLower = -0.1,
    riskDiffUpper = 0.1,
    pi1 = 0.12,
    pi2 = 0.12,
    kMax = 3,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designRiskDiffEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8007)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 369)
  testthat::expect_equal(round(design1$overallResults$information, 4), 873.5795)
  testthat::expect_equal(
    round(design1$byStageResults$efficacyBounds, 3),
    c(3.200, 2.141, 1.695)
  )
})


testthat::test_that("getDesignRepeatedANOVAContrast: Rd example numeric regression", {
  design1 <- getDesignRepeatedANOVAContrast(
    beta = 0.1,
    ngroups = 4,
    means = c(1.5, 2.5, 2, 0),
    stDev = 5,
    corr = 0.2,
    contrast = c(1, 1, 1, -3) / 3,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designRepeatedANOVAContrast")
  testthat::expect_equal(round(design1$power, 4), 0.9012)
  testthat::expect_equal(design1$n, 71)
  testthat::expect_equal(round(design1$effectsize, 3), 0.150)
  testthat::expect_equal(design1$meanContrast, 2)
})


testthat::test_that("getDesignMeanDiffCarryover: Rd example numeric regression", {
  design1 <- getDesignMeanDiffCarryover(
    beta = 0.2,
    n = NA,
    meanDiff = 0.5,
    stDev = 1,
    design = matrix(c(1, 4, 2, 3,
                      2, 1, 3, 4,
                      3, 2, 4, 1,
                      4, 3, 1, 2),
                    4, 4, byrow = TRUE),
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designMeanDiffCarryover")
  testthat::expect_equal(round(design1$power, 4), 0.8015)
  testthat::expect_equal(design1$numberOfSubjects, 70)
})


testthat::test_that("getDesignMeanDiffCarryoverEquiv: Rd example numeric regression", {
  design1 <- getDesignMeanDiffCarryoverEquiv(
    beta = 0.2,
    n = NA,
    meanDiffLower = -1.3,
    meanDiffUpper = 1.3,
    meanDiff = 0,
    stDev = 2.2,
    design = matrix(c(1, 4, 2, 3,
                      2, 1, 3, 4,
                      3, 2, 4, 1,
                      4, 3, 1, 2),
                    4, 4, byrow = TRUE),
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designMeanDiffCarryoverEquiv")
  testthat::expect_equal(round(design1$power, 4), 0.8011)
  testthat::expect_equal(design1$numberOfSubjects, 67)
})


testthat::test_that("getDesignMeanDiffXOEquiv: Rd examples numeric regression", {
  design1 <- getDesignMeanDiffXOEquiv(
    beta = 0.1,
    n = NA,
    meanDiffLower = -1.3,
    meanDiffUpper = 1.3,
    meanDiff = 0,
    stDev = 2.2,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignMeanDiffXOEquiv(
    beta = 0.1,
    n = NA,
    meanDiffLower = -1.3,
    meanDiffUpper = 1.3,
    meanDiff = 0,
    stDev = 2.2,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designMeanDiffXOEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9024)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 64)

  testthat::expect_s3_class(design2, "designMeanDiffXOEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9033)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 64)
})


testthat::test_that("getDesignMeanDiffXO: Rd example numeric regression", {
  design1 <- getDesignMeanDiffXO(
    beta = 0.2,
    n = NA,
    meanDiff = 75,
    stDev = 150,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designMeanDiffXO")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8009)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 51)
})


testthat::test_that("getDesignMeanRatioXO: Rd example numeric regression", {
  design1 <- getDesignMeanRatioXO(
    beta = 0.1,
    n = NA,
    meanRatio = 1.25,
    CV = 0.25,
    alpha = 0.05,
    normalApproximation = FALSE
  )

  testthat::expect_s3_class(design1, "designMeanRatioXO")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9078)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 23)
})


testthat::test_that("getDesignMeanRatioXOEquiv: Rd examples numeric regression", {
  design1 <- getDesignMeanRatioXOEquiv(
    beta = 0.1,
    n = NA,
    meanRatioLower = 0.8,
    meanRatioUpper = 1.25,
    meanRatio = 1,
    CV = 0.35,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignMeanRatioXOEquiv(
    beta = 0.1,
    n = NA,
    meanRatioLower = 0.8,
    meanRatioUpper = 1.25,
    meanRatio = 1,
    CV = 0.35,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designMeanRatioXOEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9033)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 52)

  testthat::expect_s3_class(design2, "designMeanRatioXOEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9024)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 52)
})


testthat::test_that("getDesignOrderedBinom: Rd example numeric regression", {
  design1 <- getDesignOrderedBinom(
    beta = 0.1,
    ngroups = 3,
    pi = c(0.1, 0.25, 0.5),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designOrderedBinom")
  testthat::expect_equal(round(design1$power, 4), 0.9011)
  testthat::expect_equal(design1$n, 75)
})


testthat::test_that("getDesignPairedMeanDiff: Rd examples numeric regression", {
  design1 <- getDesignPairedMeanDiff(
    beta = 0.1,
    n = NA,
    pairedDiffH0 = 0,
    pairedDiff = -2,
    stDev = 5,
    kMax = 5,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignPairedMeanDiff(
    beta = 0.1,
    n = NA,
    pairedDiffH0 = 0,
    pairedDiff = -2,
    stDev = 5,
    normalApproximation = FALSE,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designPairedMeanDiff")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9033)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 56)

  testthat::expect_s3_class(design2, "designPairedMeanDiff")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9016)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 68)
})


testthat::test_that("getDesignPairedMeanDiffEquiv: Rd examples numeric regression", {
  design1 <- getDesignPairedMeanDiffEquiv(
    beta = 0.1,
    n = NA,
    pairedDiffLower = -1.3,
    pairedDiffUpper = 1.3,
    pairedDiff = 0,
    stDev = 2.2,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignPairedMeanDiffEquiv(
    beta = 0.1,
    n = NA,
    pairedDiffLower = -1.3,
    pairedDiffUpper = 1.3,
    pairedDiff = 0,
    stDev = 2.2,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designPairedMeanDiffEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9024)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 32)

  testthat::expect_s3_class(design2, "designPairedMeanDiffEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9064)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 33)
})


testthat::test_that("getDesignMeanDiffMMRM: Rd example numeric regression", {
  design1 <- getDesignMeanDiffMMRM(
    beta = 0.1,
    meanDiff = 0.5,
    k = 2,
    t = c(1, 2),
    accrualIntensity = 40,
    accrualDuration = 1
  )

  testthat::expect_s3_class(design1, "designMeanDiffMMRM")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9015)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 169)
})


testthat::test_that("getDesignOneSlope: Rd example numeric regression", {
  design1 <- getDesignOneSlope(
    beta = 0.1,
    n = NA,
    slope = 1,
    stDev = 1,
    stDevCovariate = 1,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designOneSlope")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9126)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 11)
})


testthat::test_that("getDesignPairedMeanRatio: Rd example numeric regression", {
  design1 <- getDesignPairedMeanRatio(
    beta = 0.1,
    n = NA,
    pairedRatio = 1.2,
    CV = 0.35,
    alpha = 0.05,
    normalApproximation = FALSE
  )

  testthat::expect_s3_class(design1, "designPairedMeanRatio")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9069)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 32)
})


testthat::test_that("getDesignPairedMeanRatioEquiv: Rd examples numeric regression", {
  design1 <- getDesignPairedMeanRatioEquiv(
    beta = 0.1,
    n = NA,
    pairedRatioLower = 0.8,
    pairedRatioUpper = 1.25,
    pairedRatio = 1,
    CV = 0.35,
    kMax = 4,
    alpha = 0.05,
    typeAlphaSpending = "sfOF"
  )

  design2 <- getDesignPairedMeanRatioEquiv(
    beta = 0.1,
    n = NA,
    pairedRatioLower = 0.8,
    pairedRatioUpper = 1.25,
    pairedRatio = 1,
    CV = 0.35,
    normalApproximation = FALSE,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designPairedMeanRatioEquiv")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9029)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 26)

  testthat::expect_s3_class(design2, "designPairedMeanRatioEquiv")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9061)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 27)
})


testthat::test_that("getDesignRiskDiffExact: Rd example numeric regression", {
  design1 <- getDesignRiskDiffExact(
    n = 50,
    pi1 = 0.6,
    pi2 = 0.25,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(design1$power, 0.694611217040014)
  testthat::expect_equal(design1$n, 50)
})


testthat::test_that("getDesignRiskDiffExactEquiv: Rd example numeric regression", {
  design1 <- getDesignRiskDiffExactEquiv(
    n = 200,
    riskDiffLower = -0.2,
    riskDiffUpper = 0.2,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(round(design1$power, 4), 0.9147)
  testthat::expect_equal(design1$n, 200)
})


testthat::test_that("getDesignRiskRatioExact: Rd example numeric regression", {
  design1 <- getDesignRiskRatioExact(
    beta = 0.2,
    riskRatioH0 = 0.7,
    pi1 = 0.95,
    pi2 = 0.95,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(round(design1$power, 4), 0.8573)
  testthat::expect_equal(design1$n, 39)
})


testthat::test_that("getDesignRiskRatioExactEquiv: Rd example numeric regression", {
  design1 <- getDesignRiskRatioExactEquiv(
    n = 200,
    riskRatioLower = 0.8,
    riskRatioUpper = 1.25,
    pi1 = 0.775,
    pi2 = 0.775,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "data.frame")
  testthat::expect_equal(round(design1$power, 4), 0.7514)
  testthat::expect_equal(design1$n, 200)
})


testthat::test_that("getDesignRiskRatioFM: Rd example numeric regression", {
  design1 <- getDesignRiskRatioFM(
    beta = 0.2,
    riskRatioH0 = 1.3,
    pi1 = 0.125,
    pi2 = 0.125,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designRiskRatioFM")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8001)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 2531)
})


testthat::test_that("getDesignPairedPropMcNemar: Rd examples numeric regression", {
  design1 <- getDesignPairedPropMcNemar(
    beta = 0.1,
    n = NA,
    pDiscordant = 0.16,
    riskDiff = 0.1,
    alpha = 0.025
  )

  design2 <- getDesignPairedPropMcNemar(
    beta = 0.1,
    n = NA,
    pDiscordant = 0.16,
    riskDiff = 0.1,
    alpha = 0.025,
    kMax = 3,
    typeAlphaSpending = "sfOF"
  )

  testthat::expect_s3_class(design1, "designPairedPropMcNemar")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9001)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 164)

  testthat::expect_s3_class(design2, "designPairedPropMcNemar")
  testthat::expect_equal(round(design2$overallResults$overallReject, 4), 0.9016)
  testthat::expect_equal(design2$overallResults$numberOfSubjects, 167)
  testthat::expect_equal(
    round(design2$byStageResults$efficacyBounds, 3),
    c(3.698, 2.516, 1.993)
  )
})


testthat::test_that("getDesignSlopeDiff: Rd example numeric regression", {
  design1 <- getDesignSlopeDiff(
    beta = 0.1,
    n = NA,
    slopeDiff = -0.5,
    stDev = 10,
    stDevCovariate = 6,
    normalApproximation = FALSE,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designSlopeDiff")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.9000)
  testthat::expect_equal(design1$overallResults$numberOfSubjects, 469)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 1.965)
})


testthat::test_that("getDesignSlopeDiffMMRM: Rd example numeric regression", {
  design1 <- getDesignSlopeDiffMMRM(
    beta = 0.2,
    slopeDiff = log(1.15) / 52,
    stDev = sqrt(.182),
    stDevIntercept = sqrt(.238960),
    stDevSlope = sqrt(.000057),
    corrInterceptSlope = .003688 / sqrt(.238960 * .000057),
    w = 8,
    N = 10000,
    accrualIntensity = 15,
    gamma1 = 1 / (4.48 * 52),
    gamma2 = 1 / (4.48 * 52),
    accrualDuration = NA,
    followupTime = 8,
    alpha = 0.025
  )

  testthat::expect_s3_class(design1, "designSlopeDiffMMRM")
  testthat::expect_equal(round(design1$overallResults$overallReject, 4), 0.8004)
  testthat::expect_equal(round(design1$overallResults$information, 0), 1087517)
  testthat::expect_equal(round(design1$overallResults$numberOfSubjects, 0), 1013)
  testthat::expect_equal(round(design1$overallResults$studyDuration, 4), 75.5333)
  testthat::expect_equal(round(design1$overallResults$accrualDuration, 4), 67.5333)
  testthat::expect_equal(round(design1$byStageResults$efficacyBounds[1], 3), 1.960)
})


testthat::test_that("getDesignTwoMultinom: Rd example numeric regression", {
  design1 <- getDesignTwoMultinom(
    beta = 0.1,
    ncats = 3,
    pi1 = c(0.3, 0.35),
    pi2 = c(0.2, 0.3),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designTwoMultinom")
  testthat::expect_equal(round(design1$power, 4), 0.9000)
  testthat::expect_equal(design1$n, 503)
  testthat::expect_equal(round(design1$effectsize, 4), 0.0252)
})


testthat::test_that("getDesignTwoWayANOVA: Rd example numeric regression", {
  design1 <- getDesignTwoWayANOVA(
    beta = 0.1,
    nlevelsA = 2,
    nlevelsB = 2,
    means = matrix(c(0.5, 4.7, 0.4, 6.9), 2, 2, byrow = TRUE),
    stDev = 2,
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designTwoWayANOVA")
  testthat::expect_equal(design1$powerdf$n[1], 156)
  testthat::expect_equal(round(design1$powerdf$powerA[1], 4), 0.9028)
  testthat::expect_equal(round(design1$powerdf$powerB[1], 4), 1.0000)
  testthat::expect_equal(round(design1$powerdf$powerAB[1], 4), 0.9461)
  testthat::expect_equal(round(design1$effectsizeAB, 4), 0.0827)
})


testthat::test_that("getDesignUnorderedBinom: Rd example numeric regression", {
  design1 <- getDesignUnorderedBinom(
    beta = 0.1,
    ngroups = 3,
    pi = c(0.1, 0.25, 0.5),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designUnorderedBinom")
  testthat::expect_equal(round(design1$power, 4), 0.9020)
  testthat::expect_equal(design1$n, 95)
  testthat::expect_equal(round(design1$effectsize, 4), 0.1341)
})


testthat::test_that("getDesignUnorderedMultinom: Rd example numeric regression", {
  design1 <- getDesignUnorderedMultinom(
    beta = 0.1,
    ngroups = 3,
    ncats = 4,
    pi = matrix(c(0.230, 0.320, 0.272,
                  0.358, 0.442, 0.154,
                  0.142, 0.036, 0.039),
                3, 3, byrow = TRUE),
    allocationRatioPlanned = c(2, 2, 1),
    alpha = 0.05
  )

  testthat::expect_s3_class(design1, "designUnorderedMultinom")
  testthat::expect_equal(round(design1$power, 4), 0.9083)
  testthat::expect_equal(design1$n, 40)
  testthat::expect_equal(round(design1$effectsize, 4), 0.4466)
})
