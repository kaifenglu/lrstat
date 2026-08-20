testthat::test_that("set.seed reproducibility across all LR simulation entry points", {
  set.seed(101)
  sim_a <- lrsim(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    n = 132,
    plannedEvents = c(60, 120),
    maxNumberOfIterations = 80,
    seed = 101,
    nthreads = 1
  )
  set.seed(101)
  sim_b <- lrsim(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    n = 132,
    plannedEvents = c(60, 120),
    maxNumberOfIterations = 80,
    seed = 101,
    nthreads = 1
  )
  testthat::expect_equal(sim_a$overview$overallReject, sim_b$overview$overallReject)
  testthat::expect_equal(sim_a$sumdata$logRankStatistic, sim_b$sumdata$logRankStatistic)

  set.seed(102)
  multiarm_a <- lrsim_multiarm(
    M = 2,
    kMax = 3,
    criticalValues = matrix(c(3.880, 2.747, 2.275, 3.710, 2.511, 1.993), 3, 2),
    futilityBounds = c(0.074, 1.207),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.75, log(2) / 12),
    n = 700,
    plannedEvents = c(36, 72, 108),
    maxNumberOfIterations = 60,
    seed = 102,
    nthreads = 1
  )
  set.seed(102)
  multiarm_b <- lrsim_multiarm(
    M = 2,
    kMax = 3,
    criticalValues = matrix(c(3.880, 2.747, 2.275, 3.710, 2.511, 1.993), 3, 2),
    futilityBounds = c(0.074, 1.207),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.75, log(2) / 12),
    n = 700,
    plannedEvents = c(36, 72, 108),
    maxNumberOfIterations = 60,
    seed = 102,
    nthreads = 1
  )
  testthat::expect_equal(multiarm_a$overview$overallReject, multiarm_b$overview$overallReject)
  testthat::expect_equal(multiarm_a$sumdata2$logRankStatistic, multiarm_b$sumdata2$logRankStatistic)

  set.seed(103)
  seam_a <- lrsim_seamless(
    M = 2,
    K = 2,
    criticalValues = c(3.882, 2.733, 2.222),
    futilityBounds = c(0.259, 1.201),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.7, log(2) / 12),
    n = 700,
    plannedEvents = c(42, 84, 126),
    maxNumberOfIterations = 60,
    seed = 103,
    nthreads = 1
  )
  set.seed(103)
  seam_b <- lrsim_seamless(
    M = 2,
    K = 2,
    criticalValues = c(3.882, 2.733, 2.222),
    futilityBounds = c(0.259, 1.201),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.7, log(2) / 12),
    n = 700,
    plannedEvents = c(42, 84, 126),
    maxNumberOfIterations = 60,
    seed = 103,
    nthreads = 1
  )
  testthat::expect_equal(seam_a$overview$overallReject, seam_b$overview$overallReject)
  testthat::expect_equal(seam_a$sumdata2$logRankStatistic, seam_b$sumdata2$logRankStatistic)

  set.seed(104)
  sim2e_a <- lrsim2e(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30,
    n = 420,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 104,
    nthreads = 1
  )
  set.seed(104)
  sim2e_b <- lrsim2e(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30,
    n = 420,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 104,
    nthreads = 1
  )
  testthat::expect_equal(sim2e_a$sumdata$logRankStatistic, sim2e_b$sumdata$logRankStatistic)

  set.seed(105)
  sim2e3a_a <- lrsim2e3a(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12 * 0.70,
    lambda3pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30 * 0.75,
    lambda3os = log(2) / 30,
    n = 700,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 105,
    nthreads = 1
  )
  set.seed(105)
  sim2e3a_b <- lrsim2e3a(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12 * 0.70,
    lambda3pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30 * 0.75,
    lambda3os = log(2) / 30,
    n = 700,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 105,
    nthreads = 1
  )
  testthat::expect_equal(sim2e3a_a$sumdata$logRankStatistic13, sim2e3a_b$sumdata$logRankStatistic13)

  set.seed(106)
  sim3a_a <- lrsim3a(
    kMax = 3,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambda1 = log(2) / 12 * 0.60,
    lambda2 = log(2) / 12 * 0.70,
    lambda3 = log(2) / 12,
    n = 700,
    plannedEvents = c(186, 259, 295),
    maxNumberOfIterations = 50,
    seed = 106,
    nthreads = 1
  )
  set.seed(106)
  sim3a_b <- lrsim3a(
    kMax = 3,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambda1 = log(2) / 12 * 0.60,
    lambda2 = log(2) / 12 * 0.70,
    lambda3 = log(2) / 12,
    n = 700,
    plannedEvents = c(186, 259, 295),
    maxNumberOfIterations = 50,
    seed = 106,
    nthreads = 1
  )
  testthat::expect_equal(sim3a_a$sumdata$logRankStatistic13, sim3a_b$sumdata$logRankStatistic13)

  set.seed(107)
  simsub_a <- lrsimsub(
    kMax = 2,
    kMaxitt = 2,
    allocation1 = 1,
    allocation2 = 1,
    accrualTime = seq(0, 9),
    accrualIntensity = c(seq(10, 70, 10), rep(70, 3)),
    piecewiseSurvivalTime = c(0, 12, 24),
    p_pos = 0.6,
    lambda1itt = c(0.00256, 0.00383, 0.00700),
    lambda2itt = c(0.00427, 0.00638, 0.01167),
    lambda1pos = c(0.00299, 0.00430, 0.01064),
    lambda2pos = c(0.00516, 0.00741, 0.01835),
    gamma1itt = -log(1 - 0.04) / 12,
    gamma2itt = -log(1 - 0.04) / 12,
    gamma1pos = -log(1 - 0.04) / 12,
    gamma2pos = -log(1 - 0.04) / 12,
    n = 500,
    plannedEvents = c(108, 144),
    maxNumberOfIterations = 50,
    seed = 107,
    nthreads = 1
  )
  set.seed(107)
  simsub_b <- lrsimsub(
    kMax = 2,
    kMaxitt = 2,
    allocation1 = 1,
    allocation2 = 1,
    accrualTime = seq(0, 9),
    accrualIntensity = c(seq(10, 70, 10), rep(70, 3)),
    piecewiseSurvivalTime = c(0, 12, 24),
    p_pos = 0.6,
    lambda1itt = c(0.00256, 0.00383, 0.00700),
    lambda2itt = c(0.00427, 0.00638, 0.01167),
    lambda1pos = c(0.00299, 0.00430, 0.01064),
    lambda2pos = c(0.00516, 0.00741, 0.01835),
    gamma1itt = -log(1 - 0.04) / 12,
    gamma2itt = -log(1 - 0.04) / 12,
    gamma1pos = -log(1 - 0.04) / 12,
    gamma2pos = -log(1 - 0.04) / 12,
    n = 500,
    plannedEvents = c(108, 144),
    maxNumberOfIterations = 50,
    seed = 107,
    nthreads = 1
  )
  testthat::expect_equal(simsub_a$sumdata$logRankStatistic, simsub_b$sumdata$logRankStatistic)
})


testthat::test_that("simulation summaries contain expected fields and bounded probabilities", {
  sim <- lrsim(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    n = 132,
    plannedEvents = c(60, 120),
    maxNumberOfIterations = 60,
    seed = 201,
    nthreads = 1
  )
  testthat::expect_s3_class(sim, "lrsim")
  testthat::expect_true(all(c("overview", "sumdata") %in% names(sim)))
  testthat::expect_true(sim$overview$overallReject >= 0 && sim$overview$overallReject <= 1)
  testthat::expect_true(all(sim$sumdata$stageNumber %in% 1:2))

  multiarm <- lrsim_multiarm(
    M = 2,
    kMax = 3,
    criticalValues = matrix(c(3.880, 2.747, 2.275, 3.710, 2.511, 1.993), 3, 2),
    futilityBounds = c(0.074, 1.207),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.75, log(2) / 12),
    n = 700,
    plannedEvents = c(36, 72, 108),
    maxNumberOfIterations = 60,
    seed = 202,
    nthreads = 1
  )
  testthat::expect_s3_class(multiarm, "lrsim_multiarm")
  testthat::expect_true(all(c("overview", "sumdata1", "sumdata2") %in% names(multiarm)))
  testthat::expect_true(multiarm$overview$overallReject >= 0 && multiarm$overview$overallReject <= 1)
  testthat::expect_true(multiarm$overview$overallFutility >= 0 && multiarm$overview$overallFutility <= 1)

  seam <- lrsim_seamless(
    M = 2,
    K = 2,
    criticalValues = c(3.882, 2.733, 2.222),
    futilityBounds = c(0.259, 1.201),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.7, log(2) / 12),
    n = 700,
    plannedEvents = c(42, 84, 126),
    maxNumberOfIterations = 60,
    seed = 203,
    nthreads = 1
  )
  testthat::expect_s3_class(seam, "lrsim_seamless")
  testthat::expect_true(all(c("overview", "sumdata1", "sumdata2") %in% names(seam)))
  testthat::expect_true(seam$overview$overallReject >= 0 && seam$overview$overallReject <= 1)
  testthat::expect_true(seam$overview$overallFutility >= 0 && seam$overview$overallFutility <= 1)

  sim2e <- lrsim2e(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30,
    n = 420,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 204,
    nthreads = 1
  )
  testthat::expect_true(all(c("sumdata") %in% names(sim2e)))
  testthat::expect_true(all(c("endpoint", "logRankStatistic", "stageNumber") %in% names(sim2e$sumdata)))
  testthat::expect_true(setequal(sort(unique(stats::na.omit(sim2e$sumdata$endpoint))), c("PFS", "OS")))
  testthat::expect_true(all(sim2e$sumdata$stageNumber %in% 1:3))

  sim2e3a <- lrsim2e3a(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12 * 0.70,
    lambda3pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30 * 0.75,
    lambda3os = log(2) / 30,
    n = 700,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 50,
    seed = 205,
    nthreads = 1
  )
  testthat::expect_true(all(c("sumdata") %in% names(sim2e3a)))
  testthat::expect_true(all(c("endpoint", "logRankStatistic13", "logRankStatistic23", "logRankStatistic12") %in% names(sim2e3a$sumdata)))
  testthat::expect_true(setequal(sort(unique(stats::na.omit(sim2e3a$sumdata$endpoint))), c("PFS", "OS")))

  sim3a <- lrsim3a(
    kMax = 3,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambda1 = log(2) / 12 * 0.60,
    lambda2 = log(2) / 12 * 0.70,
    lambda3 = log(2) / 12,
    n = 700,
    plannedEvents = c(186, 259, 295),
    maxNumberOfIterations = 50,
    seed = 206,
    nthreads = 1
  )
  testthat::expect_true(all(c("sumdata") %in% names(sim3a)))
  testthat::expect_true(all(c("logRankStatistic13", "logRankStatistic23", "logRankStatistic12") %in% names(sim3a$sumdata)))
  testthat::expect_true(all(sim3a$sumdata$stageNumber %in% 1:3))

  simsub <- lrsimsub(
    kMax = 2,
    kMaxitt = 2,
    allocation1 = 1,
    allocation2 = 1,
    accrualTime = seq(0, 9),
    accrualIntensity = c(seq(10, 70, 10), rep(70, 3)),
    piecewiseSurvivalTime = c(0, 12, 24),
    p_pos = 0.6,
    lambda1itt = c(0.00256, 0.00383, 0.00700),
    lambda2itt = c(0.00427, 0.00638, 0.01167),
    lambda1pos = c(0.00299, 0.00430, 0.01064),
    lambda2pos = c(0.00516, 0.00741, 0.01835),
    gamma1itt = -log(1 - 0.04) / 12,
    gamma2itt = -log(1 - 0.04) / 12,
    gamma1pos = -log(1 - 0.04) / 12,
    gamma2pos = -log(1 - 0.04) / 12,
    n = 500,
    plannedEvents = c(108, 144),
    maxNumberOfIterations = 50,
    seed = 207,
    nthreads = 1
  )
  testthat::expect_true(all(c("sumdata") %in% names(simsub)))
  testthat::expect_true(all(c("population", "logRankStatistic") %in% names(simsub$sumdata)))
  testthat::expect_true(all(stats::na.omit(simsub$sumdata$population) %in% c("ITT", "Biomarker+", "Biomarker-")))
})


testthat::test_that("simulation functions reproduce seeded expected summary results", {
  sim <- lrsim(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    n = 132,
    plannedEvents = c(60, 120),
    maxNumberOfIterations = 120,
    seed = 301,
    nthreads = 1
  )
  testthat::expect_equal(round(sim$overview$overallReject, 6), 0.783333)
  testthat::expect_equal(round(sim$overview$expectedNumberOfEvents, 1), 109.5)

  multiarm <- lrsim_multiarm(
    M = 2,
    kMax = 3,
    criticalValues = matrix(c(3.880, 2.747, 2.275, 3.710, 2.511, 1.993), 3, 2),
    futilityBounds = c(0.074, 1.207),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.75, log(2) / 12),
    n = 700,
    plannedEvents = c(36, 72, 108),
    maxNumberOfIterations = 100,
    seed = 302,
    nthreads = 1
  )
  testthat::expect_equal(round(multiarm$overview$overallReject, 2), 0.89)
  testthat::expect_equal(round(multiarm$overview$overallFutility, 2), 0.11)
  testthat::expect_equal(round(unname(multiarm$overview$expectedNumberOfEvents[4]), 2), 126.05)

  seam <- lrsim_seamless(
    M = 2,
    K = 2,
    criticalValues = c(3.882, 2.733, 2.222),
    futilityBounds = c(0.259, 1.201),
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambdas = list(log(2) / 12 * 0.5, log(2) / 12 * 0.7, log(2) / 12),
    n = 700,
    plannedEvents = c(42, 84, 126),
    maxNumberOfIterations = 100,
    seed = 303,
    nthreads = 1
  )
  testthat::expect_equal(round(seam$overview$overallReject, 2), 0.87)
  testthat::expect_equal(round(seam$overview$overallFutility, 2), 0.13)
  testthat::expect_equal(round(seam$overview$expectedNumberOfEvents, 2), 114.57)

  sim2e <- lrsim2e(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30,
    n = 420,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 100,
    seed = 304,
    nthreads = 1
  )
  testthat::expect_equal(round(mean(sim2e$sumdata$logRankStatistic), 6), -3.218186)
  testthat::expect_equal(round(mean(sim2e$sumdata$totalEvents), 4), 195.2117)

  sim2e3a <- lrsim2e3a(
    kMax = 3,
    kMaxpfs = 2,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    rho_pd_os = 0,
    lambda1pfs = log(2) / 12 * 0.60,
    lambda2pfs = log(2) / 12 * 0.70,
    lambda3pfs = log(2) / 12,
    lambda1os = log(2) / 30 * 0.65,
    lambda2os = log(2) / 30 * 0.75,
    lambda3os = log(2) / 30,
    n = 700,
    plannedEvents = c(186, 259, 183),
    maxNumberOfIterations = 100,
    seed = 305,
    nthreads = 1
  )
  testthat::expect_equal(round(mean(sim2e3a$sumdata$logRankStatistic13), 6), -3.248537)
  testthat::expect_equal(round(mean(sim2e3a$sumdata$logRankStatistic23), 6), -2.219514)

  sim3a <- lrsim3a(
    kMax = 3,
    allocation1 = 2,
    allocation2 = 2,
    allocation3 = 1,
    accrualTime = c(0, 8),
    accrualIntensity = c(10, 28),
    lambda1 = log(2) / 12 * 0.60,
    lambda2 = log(2) / 12 * 0.70,
    lambda3 = log(2) / 12,
    n = 700,
    plannedEvents = c(186, 259, 295),
    maxNumberOfIterations = 100,
    seed = 306,
    nthreads = 1
  )
  testthat::expect_equal(round(mean(sim3a$sumdata$logRankStatistic13), 6), -3.954235)
  testthat::expect_equal(round(mean(sim3a$sumdata$logRankStatistic23), 5), -2.71789)

  simsub <- lrsimsub(
    kMax = 2,
    kMaxitt = 2,
    allocation1 = 1,
    allocation2 = 1,
    accrualTime = seq(0, 9),
    accrualIntensity = c(seq(10, 70, 10), rep(70, 3)),
    piecewiseSurvivalTime = c(0, 12, 24),
    p_pos = 0.6,
    lambda1itt = c(0.00256, 0.00383, 0.00700),
    lambda2itt = c(0.00427, 0.00638, 0.01167),
    lambda1pos = c(0.00299, 0.00430, 0.01064),
    lambda2pos = c(0.00516, 0.00741, 0.01835),
    gamma1itt = -log(1 - 0.04) / 12,
    gamma2itt = -log(1 - 0.04) / 12,
    gamma1pos = -log(1 - 0.04) / 12,
    gamma2pos = -log(1 - 0.04) / 12,
    n = 500,
    plannedEvents = c(108, 144),
    maxNumberOfIterations = 100,
    seed = 307,
    nthreads = 1
  )
  testthat::expect_equal(round(mean(simsub$sumdata$logRankStatistic), 6), -2.346592)
  testthat::expect_equal(round(mean(simsub$sumdata$totalEvents), 0), 84)
})


testthat::test_that("lrsim empirical rejection converges near analytic target in medium-size scenario", {
  analytic <- lrpower(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    accrualDuration = 12,
    followupTime = 35,
    fixedFollowup = FALSE
  )

  planned_time <- as.numeric(analytic$byStageResults$analysisTime)
  n_target <- as.integer(round(analytic$overallResults$expectedNumberOfSubjects))

  sim <- lrsim(
    kMax = 2,
    informationRates = c(0.5, 1),
    criticalValues = c(2.797, 1.977),
    accrualIntensity = 11,
    lambda1 = 0.018,
    lambda2 = 0.030,
    n = n_target,
    plannedTime = planned_time,
    maxNumberOfIterations = 1200,
    seed = 444,
    nthreads = 1
  )

  target <- analytic$overallResults$overallReject
  empirical <- sim$overview$overallReject

  testthat::expect_equal(n_target, 132)
  testthat::expect_equal(round(target, 6), 0.628949)
  testthat::expect_lt(abs(empirical - target), 0.08)
})

