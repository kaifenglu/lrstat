testthat::test_that("adaptDesign: workflow object contract and Rd example regression", {
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

  des2 <- adaptDesign(
    betaNew = 0.1,
    L = L,
    zL = zL,
    theta = 10,
    IMax = n / (4 * sigma^2),
    kMax = 3,
    informationRates = s1,
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    MullerSchafer = TRUE,
    kNew = 2,
    typeAlphaSpendingNew = "sfP"
  )

  testthat::expect_s3_class(des2, "adaptDesign")
  testthat::expect_named(des2, c("primaryTrial", "secondaryTrial", "integratedTrial"))

  testthat::expect_equal(des2$primaryTrial$L, 1L)
  testthat::expect_equal(des2$secondaryTrial$kMax, 2L)
  testthat::expect_equal(des2$integratedTrial$kMax, 3L)
  testthat::expect_true(des2$primaryTrial$MullerSchafer)
  testthat::expect_equal(round(des2$primaryTrial$conditionalPower, 4), 0.4957)
  testthat::expect_equal(round(des2$primaryTrial$predictivePower, 4), 0.3877)
  testthat::expect_equal(round(des2$secondaryTrial$maxInformation, 4), 0.1041)
  testthat::expect_equal(round(des2$integratedTrial$maxInformation, 4), 0.1199)
  testthat::expect_equal(
    round(des2$secondaryTrial$efficacyBounds, 4),
    c(1.9876, 2.0179)
  )
  testthat::expect_equal(
    round(des2$integratedTrial$efficacyBounds, 4),
    c(3.7060, 2.1820, 2.2121)
  )
  testthat::expect_true(all(diff(des2$integratedTrial$information) > 0))
})


testthat::test_that("adaptDesign_multiarm and adaptDesign_seamless: workflow objects contain required decision data", {
  des_multiarm <- adaptDesign_multiarm(
    betaNew = 0.2,
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

  des_seam <- adaptDesign_seamless(
    betaNew = 0.1,
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

  testthat::expect_s3_class(des_multiarm, "adaptDesign_multiarm")
  testthat::expect_named(des_multiarm, c("primaryTrial", "secondaryTrial", "integratedTrial"))
  testthat::expect_equal(des_multiarm$secondaryTrial$selected, 2L)
  testthat::expect_equal(round(des_multiarm$primaryTrial$conditionalAlpha, 4), 0.0597)
  testthat::expect_equal(round(des_multiarm$primaryTrial$conditionalPower, 4), 0.4955)
  testthat::expect_equal(round(des_multiarm$secondaryTrial$maxInformation, 4), 93.2245)
  testthat::expect_equal(round(des_multiarm$integratedTrial$maxInformation, 4), 133.7246)
  testthat::expect_equal(nrow(des_multiarm$secondaryTrial$byHypothesisBounds), 1)
  testthat::expect_equal(nrow(des_multiarm$integratedTrial$byIntersectionBounds), 2)
  testthat::expect_true(all(des_multiarm$secondaryTrial$cumulativeAlphaSpent >= 0))

  testthat::expect_s3_class(des_seam, "adaptDesign_seamless")
  testthat::expect_named(des_seam, c("primaryTrial", "secondaryTrial", "integratedTrial"))
  testthat::expect_equal(des_seam$primaryTrial$K, 2L)
  testthat::expect_equal(round(des_seam$primaryTrial$conditionalPower, 4), 0.4402)
  testthat::expect_equal(round(des_seam$secondaryTrial$maxInformation, 4), 49.5108)
  testthat::expect_equal(round(des_seam$integratedTrial$maxInformation, 4), 69.5108)
  testthat::expect_equal(
    round(des_seam$integratedTrial$efficacyBounds, 4),
    c(3.8521, 2.7238, 2.0741)
  )
  testthat::expect_true(all(diff(des_seam$integratedTrial$information) > 0))
})


testthat::test_that("exit probability functions: probabilities are bounded and terminal outcomes sum correctly", {
  out <- exitprob(
    b = c(2.963, 2.359, 2.014),
    a = c(-0.264, 0.599, 2.014),
    theta = c(0.141, 0.204, 0.289),
    I = c(81, 121, 160)
  )

  I_multiarm <- 95 / 2 * seq(1, 3) / 3
  b_multiarm <- c(3.886562, 2.748214, 2.243907)
  p0_multiarm <- exitprob_multiarm(M = 2, theta = c(0, 0), kMax = 3, b = b_multiarm, I = I_multiarm)
  p1_multiarm <- exitprob_multiarm(M = 2, theta = c(0.3, 0.5), kMax = 3, b = b_multiarm, I = I_multiarm)

  I_seam <- 110 / 2 * seq(1, 3) / 3
  b_seam <- c(3.776605, 2.670463, 2.180424)
  a_seam <- c(0, 0.5, b_seam[3])
  p0_seam <- exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = b_seam, I = I_seam)
  p1_seam <- exitprob_seamless(M = 2, theta = c(0.3, 0.5), K = 2, b = b_seam, a = a_seam, I = I_seam)

  testthat::expect_named(out, c("exitProbUpper", "exitProbLower"))
  testthat::expect_equal(
    round(out$exitProbUpper, 4),
    c(0.0451, 0.4091, 0.4459)
  )
  testthat::expect_equal(
    round(out$exitProbLower, 4),
    c(0.0626, 0.0207, 0.0165)
  )
  testthat::expect_true(all(out$exitProbUpper >= 0 & out$exitProbUpper <= 1))
  testthat::expect_true(all(out$exitProbLower >= 0 & out$exitProbLower <= 1))
  testthat::expect_equal(sum(out$exitProbUpper) + sum(out$exitProbLower), 1, tolerance = 1e-8)

  testthat::expect_named(p0_multiarm, c("exitProbUpper", "exitProbLower"))
  testthat::expect_true(all(p0_multiarm$exitProbUpper >= 0 & p0_multiarm$exitProbUpper <= 1))
  testthat::expect_true(all(p1_multiarm$exitProbUpper >= 0 & p1_multiarm$exitProbUpper <= 1))
  testthat::expect_true(sum(p1_multiarm$exitProbUpper) > sum(p0_multiarm$exitProbUpper))
  testthat::expect_true(all(cumsum(p1_multiarm$exitProbUpper) <= 1 + 1e-10))

  testthat::expect_named(
    p1_seam,
    c("exitProbUpper", "exitProbLower", "exitProbByArmUpper",
      "exitProbByArmLower", "selectionProb")
  )
  testthat::expect_equal(dim(p0_seam$exitProbByArmUpper), c(3, 2))
  testthat::expect_equal(dim(p1_seam$exitProbByArmLower), c(3, 2))
  testthat::expect_true(all(p1_seam$exitProbUpper >= 0 & p1_seam$exitProbUpper <= 1))
  testthat::expect_true(all(p1_seam$exitProbLower >= 0 & p1_seam$exitProbLower <= 1))
  testthat::expect_equal(sum(p1_seam$selectionProb), 1, tolerance = 1e-8)
  testthat::expect_equal(
    sum(p1_seam$exitProbUpper) + sum(p1_seam$exitProbLower),
    1,
    tolerance = 1e-6
  )
})


testthat::test_that("adaptive workflow and exit probability functions: invalid stage rules are rejected", {
  testthat::expect_error(
    adaptDesign(
      betaNew = 0.1,
      L = 3,
      zL = 1.5,
      theta = 1,
      IMax = 10,
      kMax = 3,
      informationRates = c(1 / 3, 2 / 3, 1),
      alpha = 0.025
    ),
    regexp = "kMax must be greater than L|greater than L"
  )

  testthat::expect_error(
    adaptDesign_multiarm(
      betaNew = 0.2,
      M = 2,
      r = 1,
      corr_known = FALSE,
      L = 1,
      zL = c(1.2, 1.3),
      theta = c(0.1, 0.2),
      IMax = 10,
      kMax = 2,
      informationRates = c(1 / 2, 1),
      alpha = 0.025,
      MNew = 1,
      selected = 3,
      rNew = 1
    ),
    regexp = "Invalid value in selected|selected"
  )

  testthat::expect_error(
    adaptDesign_seamless(
      betaNew = 0.1,
      M = 2,
      r = 1,
      corr_known = FALSE,
      L = 2,
      zL = 1.5,
      theta = 0.2,
      IMax = 10,
      K = 2,
      informationRates = c(1 / 3, 2 / 3, 1),
      alpha = 0.025,
      kNew = 1
    ),
    regexp = "K must be greater than L|greater than L"
  )

  testthat::expect_error(
    exitprob(b = c(2, 1.5), I = c(2, 1)),
    regexp = "increasing"
  )

  testthat::expect_error(
    exitprob_multiarm(M = 2, theta = c(0, 0), kMax = 3, b = c(3, 2, 1), I = c(1, 1, 2)),
    regexp = "increasing"
  )

  testthat::expect_error(
    exitprob_seamless(M = 2, theta = c(0, 0), K = 2, b = c(3, 2, 1), I = c(1, 1, 2)),
    regexp = "increasing"
  )
})
