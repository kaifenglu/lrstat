testthat::test_that("kmpower and kmpower1s: power increases with sample size and effect size", {
  km <- kmpower(
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

  km_big <- kmpower(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 2 * (26 / 9 * seq(1, 9)),
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

  km_weaker <- kmpower(
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
    lambda1 = c(0.0533, 0.0450, 1.5 * 0.0533, 1.5 * 0.0450),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  km1s <- kmpower1s(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survH0 = 0.30,
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

  km1s_big <- kmpower1s(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survH0 = 0.30,
    accrualTime = seq(0, 8),
    accrualIntensity = 2 * (26 / 9 * seq(1, 9)),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    gamma = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )

  testthat::expect_s3_class(km, "kmpower")
  testthat::expect_named(km, c("byStageResults", "overallResults", "settings", "byTreatmentCounts"))
  testthat::expect_equal(round(km$overallResults$overallReject, 4), 0.7634)
  testthat::expect_equal(km$overallResults$numberOfSubjects, 468)
  testthat::expect_equal(km$overallResults$studyDuration, 40)

  testthat::expect_true(km_big$overallResults$overallReject > km$overallResults$overallReject)
  testthat::expect_true(km$overallResults$overallReject > km_weaker$overallResults$overallReject)
  testthat::expect_true(km1s_big$overallResults$overallReject > km1s$overallResults$overallReject)

  testthat::expect_s3_class(km1s, "kmpower1s")
  testthat::expect_named(km1s, c("byStageResults", "overallResults", "settings"))
  testthat::expect_equal(round(km1s$overallResults$overallReject, 4), 0.9557)
})


testthat::test_that("kmsamplesize family: solved designs invert target power within tolerance", {
  ss <- kmsamplesize(
    beta = 0.25,
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
    followupTime = NA,
    fixedFollowup = FALSE
  )

  ss1s <- kmsamplesize1s(
    beta = 0.2,
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survH0 = 0.30,
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

  ss_eq <- kmsamplesizeequiv(
    beta = 0.1,
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survDiffLower = -0.13,
    survDiffUpper = 0.13,
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

  testthat::expect_named(ss, c("resultsUnderH1", "resultsUnderH0"))
  testthat::expect_s3_class(ss$resultsUnderH1, "kmpower")
  testthat::expect_equal(round(ss$resultsUnderH1$overallResults$overallReject, 4), 0.7500)
  testthat::expect_equal(round(ss$resultsUnderH0$overallResults$overallReject, 4), 0.0250)
  testthat::expect_equal(round(ss$resultsUnderH1$overallResults$followupTime, 4), 14.4496)
  testthat::expect_equal(ss$resultsUnderH1$overallResults$numberOfSubjects, 468)

  testthat::expect_named(ss1s, c("resultsUnderH1", "resultsUnderH0"))
  testthat::expect_s3_class(ss1s$resultsUnderH1, "kmpower1s")
  testthat::expect_equal(round(ss1s$resultsUnderH1$overallResults$overallReject, 4), 0.8000)
  testthat::expect_equal(round(ss1s$resultsUnderH1$overallResults$followupTime, 4), 4.5323)
  testthat::expect_equal(ss1s$resultsUnderH1$overallResults$numberOfSubjects, 468)

  testthat::expect_s3_class(ss_eq, "kmpowerequiv")
  testthat::expect_equal(round(ss_eq$overallResults$overallReject, 4), 0.9000)
  testthat::expect_equal(round(ss_eq$overallResults$accrualDuration, 4), 23.9615)
  testthat::expect_equal(round(ss_eq$overallResults$followupTime, 4), 17.6773)
  testthat::expect_equal(ss_eq$overallResults$numberOfSubjects, 519)
})


testthat::test_that("one-sided and equivalence KM functions show expected directionality", {
  km_eq <- kmpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survDiffLower = -0.13,
    survDiffUpper = 0.13,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
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

  km_eq_off <- kmpowerequiv(
    kMax = 2,
    informationRates = c(0.5, 1),
    alpha = 0.05,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survDiffLower = -0.13,
    survDiffUpper = 0.13,
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

  km_two_sided <- kmpower(
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

  km_one_sided <- kmpower1s(
    kMax = 2,
    informationRates = c(0.8, 1),
    alpha = 0.025,
    typeAlphaSpending = "sfOF",
    milestone = 18,
    survH0 = 0.30,
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

  testthat::expect_s3_class(km_eq, "kmpowerequiv")
  testthat::expect_equal(round(km_eq$overallResults$overallReject, 4), 0.8609)
  testthat::expect_true(km_eq$overallResults$overallReject > km_eq_off$overallResults$overallReject)
  testthat::expect_true(km_one_sided$overallResults$overallReject > km_two_sided$overallResults$overallReject)
})



testthat::test_that("kmest: estimate and standard error", {
  library(survival)
  df1 <- kmest(aml, stratum="x", time="time", event="status",
               conftype="none")

  df2 <- summary(survfit(Surv(time, status) ~ x, data=aml,
                         conf.type="none"))

  testthat::expect_equal(df1$surv[df1$time > 0], df2$surv)
  testthat::expect_equal(df1$sesurv[df1$time > 0], df2$std.err)
})


testthat::test_that("kmest: confidence interval", {
  library(survival)
  for (conftype in c("plain", "log", "log-log", "arcsin")) {
    df1 <- kmest(aml, stratum="x", time="time", event="status",
                 conftype=conftype)

    df2 <- summary(survfit(Surv(time, status) ~ x, data=aml,
                           conf.type=conftype))

    testthat::expect_equal(df1$lower[df1$time > 0], df2$lower)
    testthat::expect_equal(df1$upper[df1$time > 0], df2$upper)
  }
})



testthat::test_that("rmdiff: unstratified test of rmst difference", {
  library(survival)
  df1 <- rmdiff(veteran, treat="trt", time="time",
                event="status", milestone=90)

  # the following values are obtained from SAS PROC LIFETEST:
  #   proc lifetest data=veteran rmst(tau=90);
  #     time time*status(0);
  #     strata trt;
  #   run;
  rmst1 = 62.75618
  rmst2 = 56.45116
  stderr1 = 4.0321
  stderr2 = 3.9790
  rmstdiffchisq = 1.2388
  pvalue = 0.2657

  testthat::expect_equal(c(round(df1$rmst1, 5),
                           round(df1$rmst2, 5),
                           round(sqrt(df1$vrmst1), 4),
                           round(sqrt(df1$vrmst2), 4),
                           round(df1$rmstDiffZ^2, 4),
                           round(df1$rmstDiffPValue, 4)),
                         c(rmst1, rmst2, stderr1, stderr2,
                           rmstdiffchisq, pvalue))
})


testthat::test_that("rmdiff: stratified test of rmst difference", {
  library(survival)
  df1 <- rmdiff(veteran, stratum="celltype", treat="trt",
                time="time", event="status", milestone=90)

  # Of note, the stratified results are different from SAS PROC LIFETEST:
  #   proc lifetest data=veteran rmst(tau=90);
  #     time time*status(0);
  #     strata celltype / group = trt;
  #   run;
  # This is because SAS adds up the rmst diffs and variances across strata,
  # while we use the number of subjects as the weight for each stratum.
  # Our approach yields a more meaningful rmst diff across strata.
  # To reproduce our results, we combine the stratum-specific information
  # from SAS PROC LIFETEST using the sample size weights across strata
  n = c(15, 20, 30, 18, 9, 18, 15, 12)
  rmst = c(68.55152, 65.95, 50.86667, 45.16667, 56.44444, 51.91667,
           84.8, 64.25)
  stderr = c(8.3317, 7.5975, 5.7978, 7.8103, 12.7986, 6.9569,
             5.0237, 8.0847)
  a = c(1, 3, 5, 7)
  b = c(2, 4, 6, 8)
  ns = n[a] + n[b]
  rmstDiffs = rmst[a] - rmst[b]
  vrmstDiffs = stderr[a]^2 + stderr[b]^2
  w = ns/sum(ns)
  rmstDiff = sum(w*rmstDiffs)
  sermstDiff = sqrt(sum(w*w*vrmstDiffs))
  rmstDiffZ = rmstDiff/sermstDiff

  testthat::expect_equal(c(round(df1$rmstDiff, 4),
                           round(df1$sermstDiff, 3),
                           round(df1$rmstDiffZ, 3)),
                         c(round(rmstDiff, 4),
                           round(sermstDiff, 3),
                           round(rmstDiffZ, 3)))
})
