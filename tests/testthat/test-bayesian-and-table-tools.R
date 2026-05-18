testthat::test_that("table outputs have expected dimensions labels and boundary behavior", {
  boin <- BOINTable(nMax = 18, pT = 0.3, phi1 = 0.6 * 0.3, phi2 = 1.4 * 0.3)

  testthat::expect_s3_class(boin, "BOINTable")
  testthat::expect_named(boin, c("settings", "decisionDataFrame", "decisionMatrix"))
  testthat::expect_equal(dim(boin$decisionDataFrame), c(189L, 3L))
  testthat::expect_equal(dim(boin$decisionMatrix), c(19L, 18L))
  testthat::expect_equal(colnames(boin$decisionDataFrame), c("n", "y", "decision"))
  testthat::expect_equal(boin$settings$lambda1, 0.23649071181459, tolerance = 1e-6)
  testthat::expect_equal(boin$settings$lambda2, 0.358519502640613, tolerance = 1e-6)
  testthat::expect_equal(boin$decisionDataFrame$decision[boin$decisionDataFrame$n == 3 & boin$decisionDataFrame$y == 1], "S")
  testthat::expect_equal(boin$decisionDataFrame$decision[boin$decisionDataFrame$n == 3 & boin$decisionDataFrame$y == 3], "DU")

  mtpi <- mTPI2Table(nMax = 18, pT = 0.3, epsilon1 = 0.05, epsilon2 = 0.05)

  testthat::expect_s3_class(mtpi, "mTPI2Table")
  testthat::expect_named(mtpi, c("settings", "subintervals", "decisionDataFrame", "decisionMatrix"))
  testthat::expect_equal(dim(mtpi$decisionDataFrame), c(189L, 3L))
  testthat::expect_equal(dim(mtpi$decisionMatrix), c(19L, 18L))
  testthat::expect_equal(colnames(mtpi$subintervals), c("lower", "upper", "decision"))
  testthat::expect_equal(mtpi$subintervals$decision[4], "S")
  testthat::expect_equal(mtpi$decisionDataFrame$decision[mtpi$decisionDataFrame$n == 2 & mtpi$decisionDataFrame$y == 2], "DU")

  s2 <- simon2stage(0.05, 0.2, 0.1, 0.3)
  testthat::expect_s3_class(s2, "data.frame")
  testthat::expect_true(all(c("n", "n1", "r1", "r", "EN0", "attainedAlpha", "attainedPower", "PET0", "design") %in% names(s2)))
  testthat::expect_equal(nrow(s2), 4L)
  testthat::expect_equal(s2$design, c("Minimax", "Admissible", "Admissible", "Optimal"))
  testthat::expect_true(all(s2$w_lower <= s2$w_upper))
})


testthat::test_that("Bayesian posterior and decision outputs match known mini examples", {
  s2 <- simon2stage(0.05, 0.2, 0.1, 0.3)
  testthat::expect_equal(s2$n, c(25, 26, 27, 29))
  testthat::expect_equal(s2$n1, c(15, 12, 11, 10))
  testthat::expect_equal(s2$r1, c(1, 1, 1, 1))
  testthat::expect_equal(s2$r, c(5, 5, 5, 5))
  testthat::expect_equal(s2$EN0, c(19.5095699875886, 16.7739662128507, 15.8422872114525, 15.0141220000001), tolerance = 1e-5)
  testthat::expect_equal(s2$attainedAlpha, c(0.0328086686626304, 0.0359671543394388, 0.0395005247254554, 0.0470863078239432), tolerance = 1e-5)
  testthat::expect_equal(s2$attainedPower, c(0.801700606072, 0.804780384089326, 0.806195368230603, 0.805062925728), tolerance = 1e-5)
  testthat::expect_equal(s2$PET0, c(0.549043, 0.659002251, 0.6973568802, 0.7360989291), tolerance = 1e-5)

  a <- simonBayesAnalysis(
    nstrata = 10,
    r = c(8, 0, 1, 1, 6, 2, 0, 0, 3, 3),
    n = c(19, 10, 26, 8, 14, 7, 8, 5, 4, 14),
    lambda = 0.5,
    gamma = 0.33,
    phi = 0.35,
    plo = 0.15
  )

  testthat::expect_equal(dim(a$case), c(1024L, 10L))
  testthat::expect_equal(length(a$prior_case), 1024L)
  testthat::expect_equal(length(a$post_case), 1024L)
  testthat::expect_equal(length(a$prior_stratum), 10L)
  testthat::expect_equal(length(a$post_stratum), 10L)
  testthat::expect_equal(sum(a$prior_case), 1, tolerance = 1e-12)
  testthat::expect_equal(sum(a$post_case), 1, tolerance = 1e-12)
  testthat::expect_equal(
    a$post_stratum,
    c(0.950120629975544, 0.0323899248727677, 0.00145815286061968,
      0.148349240807383, 0.895742212725438, 0.408975722590583,
      0.0540918875145562, 0.113259674668212, 0.820611914096308,
      0.244634489891305),
    tolerance = 1e-5
  )
})


testthat::test_that("simulation functions are reproducible under set.seed", {
  sim1 <- simonBayesSim(
    p = c(0.25, 0.25, 0.05),
    accrualIntensity = 5,
    stratumFraction = c(1 / 3, 1 / 3, 1 / 3),
    lambda = 0.33,
    gamma = 0.5,
    phi = 0.25,
    plo = 0.05,
    T = 0.8,
    maxSubjects = 50,
    plannedSubjects = seq(5, 50, 5),
    maxNumberOfIterations = 120,
    maxNumberOfRawDatasets = 1,
    seed = 314159
  )

  sim2 <- simonBayesSim(
    p = c(0.25, 0.25, 0.05),
    accrualIntensity = 5,
    stratumFraction = c(1 / 3, 1 / 3, 1 / 3),
    lambda = 0.33,
    gamma = 0.5,
    phi = 0.25,
    plo = 0.05,
    T = 0.8,
    maxSubjects = 50,
    plannedSubjects = seq(5, 50, 5),
    maxNumberOfIterations = 120,
    maxNumberOfRawDatasets = 1,
    seed = 314159
  )

  testthat::expect_named(sim1, c("rawdata", "sumdata1", "sumdata2", "overview"))
  testthat::expect_equal(dim(sim1$rawdata), c(75L, 6L))
  testthat::expect_equal(dim(sim1$sumdata1), c(1902L, 10L))
  testthat::expect_equal(dim(sim1$sumdata2), c(120L, 9L))

  testthat::expect_equal(sim1$overview$numberOfStrata, 3)
  testthat::expect_equal(sim1$overview$n_active_strata, 2)
  testthat::expect_equal(sim1$overview$true_positive, 1.59166666666667, tolerance = 1e-10)
  testthat::expect_equal(sim1$overview$false_negative, 0.4, tolerance = 1e-12)
  testthat::expect_equal(sim1$overview$false_positive, 0.15, tolerance = 1e-12)
  testthat::expect_equal(sim1$overview$true_negative, 0.833333333333333, tolerance = 1e-10)
  testthat::expect_equal(sim1$overview$n_indet_strata, 0.025, tolerance = 1e-12)
  testthat::expect_equal(sim1$overview$numberOfSubjects, 26.4166666666667, tolerance = 1e-10)

  testthat::expect_equal(sim1$sumdata2, sim2$sumdata2)
  testthat::expect_equal(sim1$overview, sim2$overview)
  testthat::expect_equal(sim1$rawdata, sim2$rawdata)

  sim3 <- simonBayesSim(
    p = c(0.25, 0.25, 0.05),
    accrualIntensity = 5,
    stratumFraction = c(1 / 3, 1 / 3, 1 / 3),
    lambda = 0.33,
    gamma = 0.5,
    phi = 0.25,
    plo = 0.05,
    T = 0.8,
    maxSubjects = 50,
    plannedSubjects = seq(5, 50, 5),
    maxNumberOfIterations = 120,
    maxNumberOfRawDatasets = 1,
    seed = 314160
  )

  testthat::expect_false(isTRUE(all.equal(sim1$sumdata2, sim3$sumdata2, tolerance = 0)))
})
