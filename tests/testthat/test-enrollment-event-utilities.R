testthat::test_that("accrual: numeric correctness for uniform and piecewise inputs", {
  # Uniform accrual from Rd example.
  uniform_t3 <- accrual(time = 3, accrualTime = 0, accrualIntensity = 20,
                        accrualDuration = 12)
  uniform_t12 <- accrual(time = 12, accrualTime = 0, accrualIntensity = 20,
                         accrualDuration = 12)

  testthat::expect_equal(uniform_t3, 60)
  testthat::expect_equal(uniform_t12, 240)

  # Piecewise accrual from Rd example.
  piecewise_vec <- accrual(time = c(2, 9), accrualTime = c(0, 3),
                           accrualIntensity = c(10, 20), accrualDuration = 12)
  piecewise_t3 <- accrual(time = 3, accrualTime = c(0, 3),
                          accrualIntensity = c(10, 20), accrualDuration = 12)

  testthat::expect_equal(piecewise_vec, c(20, 150), tolerance = 1e-10)
  testthat::expect_equal(piecewise_t3, 30, tolerance = 1e-10)

  # Three-interval regression check with hand-computed expected value.
  piecewise_t10 <- accrual(time = 10, accrualTime = c(0, 3, 7),
                           accrualIntensity = c(5, 10, 20), accrualDuration = 15)

  testthat::expect_equal(piecewise_t10, 115, tolerance = 1e-10)
})


testthat::test_that("accrual: time is clamped to valid enrollment window", {
  at_zero <- accrual(time = 0, accrualTime = 0, accrualIntensity = 20,
                     accrualDuration = 12)
  below_zero <- accrual(time = -5, accrualTime = 0, accrualIntensity = 20,
                        accrualDuration = 12)
  at_duration <- accrual(time = 12, accrualTime = 0, accrualIntensity = 20,
                         accrualDuration = 12)
  beyond_duration <- accrual(time = 20, accrualTime = 0, accrualIntensity = 20,
                             accrualDuration = 12)

  testthat::expect_equal(at_zero, 0)
  testthat::expect_equal(below_zero, 0)
  testthat::expect_equal(at_duration, 240)
  testthat::expect_equal(beyond_duration, at_duration)
})


testthat::test_that("accrual: output structure contract", {
  times <- c(1, 5, 10)
  out_vec <- accrual(time = times, accrualTime = c(0, 3),
                     accrualIntensity = c(10, 20), accrualDuration = 12)
  out_scalar <- accrual(time = 3, accrualTime = c(0, 3),
                        accrualIntensity = c(10, 20), accrualDuration = 12)

  testthat::expect_true(is.numeric(out_vec))
  testthat::expect_length(out_vec, length(times))
  testthat::expect_length(out_scalar, 1)
  testthat::expect_true(all(out_vec >= 0))

  # Non-decreasing expected accrual in calendar time.
  testthat::expect_true(all(diff(out_vec) >= 0))
})


testthat::test_that("getAccrualDurationFromN: reproduces documented piecewise example", {
  out <- getAccrualDurationFromN(
    nsubjects = c(20, 150),
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20)
  )

  testthat::expect_equal(out, c(2, 9), tolerance = 1e-10)
})


testthat::test_that("getAccrualDurationFromN: inverts accrual over piecewise schedule", {
  n <- c(5, 20, 30, 55, 95, 150)
  t <- getAccrualDurationFromN(
    nsubjects = n,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20)
  )

  n_back <- accrual(
    time = t,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    accrualDuration = max(t)
  )

  testthat::expect_equal(n_back, n, tolerance = 1e-10)
})


testthat::test_that("getAccrualDurationFromN: structure and monotonicity contract", {
  n <- c(1, 30, 31, 100)
  out <- getAccrualDurationFromN(
    nsubjects = n,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20)
  )

  testthat::expect_true(is.numeric(out))
  testthat::expect_length(out, length(n))
  testthat::expect_equal(out[2], 3, tolerance = 1e-10)
  testthat::expect_true(all(diff(out) >= 0))

  # Scalar input should return length-1 output.
  out_scalar <- getAccrualDurationFromN(
    nsubjects = 20,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20)
  )

  testthat::expect_length(out_scalar, 1)
  testthat::expect_equal(out_scalar, 2, tolerance = 1e-10)
})


testthat::test_that("patrisk and pevent: match exponential closed-form", {
  times <- c(0, 3, 9)
  lambda <- 0.0533
  gamma <- -log(1 - 0.05) / 12
  theta <- lambda + gamma

  p_risk <- patrisk(
    time = times,
    piecewiseSurvivalTime = 0,
    lambda = lambda,
    gamma = gamma
  )
  p_event <- pevent(
    time = times,
    piecewiseSurvivalTime = 0,
    lambda = lambda,
    gamma = gamma
  )

  testthat::expect_equal(p_risk, exp(-theta * times), tolerance = 1e-10)
  testthat::expect_equal(
    p_event,
    lambda / theta * (1 - exp(-theta * times)),
    tolerance = 1e-10
  )
  testthat::expect_true(all(p_event + p_risk <= 1 + 1e-12))
})


testthat::test_that("patrisk and pevent: probability bounds and monotonicity", {
  times <- seq(0, 24, by = 3)

  p_risk <- patrisk(
    time = times,
    piecewiseSurvivalTime = c(0, 6),
    lambda = c(0.0533, 0.0309),
    gamma = -log(1 - 0.05) / 12
  )
  p_event <- pevent(
    time = times,
    piecewiseSurvivalTime = c(0, 6),
    lambda = c(0.0533, 0.0309),
    gamma = -log(1 - 0.05) / 12
  )

  testthat::expect_true(all(p_risk >= 0 & p_risk <= 1))
  testthat::expect_true(all(p_event >= 0 & p_event <= 1))
  testthat::expect_true(all(diff(p_risk) <= 1e-12))
  testthat::expect_true(all(diff(p_event) >= -1e-12))
})


testthat::test_that("natrisk: agrees with allocation-accrual-probability identity", {
  t <- c(1, 5, 9, 12)
  analysis_time <- 15

  out <- natrisk(
    t = t,
    allocationRatioPlanned = 1,
    accrualTime = 0,
    accrualIntensity = 20,
    piecewiseSurvivalTime = 0,
    lambda1 = 0.1,
    lambda2 = 0.1,
    gamma1 = 0.02,
    gamma2 = 0.02,
    accrualDuration = 12,
    maxFollowupTime = 30,
    time = analysis_time
  )

  t1 <- pmin(pmin(t, 30), analysis_time)
  expected_each <- 0.5 * accrual(
    time = analysis_time - t1,
    accrualTime = 0,
    accrualIntensity = 20,
    accrualDuration = 12
  ) * patrisk(
    time = t1,
    piecewiseSurvivalTime = 0,
    lambda = 0.1,
    gamma = 0.02
  )

  testthat::expect_equal(dim(out), c(length(t), 2))
  testthat::expect_equal(out[, 1], expected_each, tolerance = 1e-10)
  testthat::expect_equal(out[, 2], expected_each, tolerance = 1e-10)
})


testthat::test_that("natrisk: returns zero when analysis time exceeds calendar time", {
  out <- natrisk(
    t = 10,
    allocationRatioPlanned = 1,
    accrualTime = 0,
    accrualIntensity = 20,
    piecewiseSurvivalTime = 0,
    lambda1 = 0.1,
    lambda2 = 0.1,
    gamma1 = 0.02,
    gamma2 = 0.02,
    accrualDuration = 12,
    maxFollowupTime = 30,
    time = 5
  )

  testthat::expect_equal(as.numeric(out), c(0, 0), tolerance = 1e-12)
})


testthat::test_that("nevent: zero hazard gives zero events", {
  times <- c(0, 6, 12, 24)

  out <- nevent(
    time = times,
    allocationRatioPlanned = 1,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0, 0),
    lambda2 = c(0, 0),
    gamma1 = 0.01,
    gamma2 = 0.01,
    accrualDuration = 12,
    maxFollowupTime = 30
  )

  testthat::expect_equal(dim(out), c(length(times), 2))
  testthat::expect_equal(out, matrix(0, nrow = length(times), ncol = 2),
                         tolerance = 1e-12)
})


testthat::test_that("nevent: monotone in calendar time and symmetric when groups match", {
  times <- seq(0, 30, by = 6)

  out <- nevent(
    time = times,
    allocationRatioPlanned = 1,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0309),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 12,
    maxFollowupTime = 30
  )

  enrolled <- accrual(
    time = times,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    accrualDuration = 12
  )

  testthat::expect_true(all(diff(out[, 1]) >= -1e-10))
  testthat::expect_true(all(diff(out[, 2]) >= -1e-10))
  testthat::expect_equal(out[, 1], out[, 2], tolerance = 1e-10)
  testthat::expect_true(all(out[, 1] <= 0.5 * enrolled + 1e-8))
  testthat::expect_true(all(out[, 2] <= 0.5 * enrolled + 1e-8))
})


# Rd example equality checks ---------------------------------------------------
# Reference values captured from the development source via pkgload::load_all().

testthat::test_that("patrisk: reproduces documented example values", {
  # Rd example: piecewise exponential survival with hazard 0.0533 in [0, 6) and
  # 0.0309 thereafter, with 5% dropout by the end of 1 year.
  out <- patrisk(
    time = c(3, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda = c(0.0533, 0.0309),
    gamma = -log(1 - 0.05) / 12
  )

  testthat::expect_equal(out, c(0.841370369971862, 0.637009970796535),
                         tolerance = 1e-10)
})


testthat::test_that("pevent: reproduces documented example values", {
  # Rd example: same piecewise exponential survival and dropout as above.
  out <- pevent(
    time = c(3, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda = c(0.0533, 0.0309),
    gamma = -log(1 - 0.05) / 12
  )

  testthat::expect_equal(out, c(0.146852650315469, 0.332689106319710),
                         tolerance = 1e-10)
})


testthat::test_that("natrisk: reproduces documented example values", {
  # Rd example: piecewise accrual, piecewise exponential survivals, 5% dropout.
  out <- natrisk(
    t = c(9, 24),
    allocationRatioPlanned = 1,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 12,
    maxFollowupTime = 30,
    time = 30
  )

  testthat::expect_equal(
    out,
    matrix(c(66.8860469336362, 16.9128874275027,
             62.5390010776441, 11.3008269108836),
           nrow = 2, ncol = 2),
    tolerance = 1e-8
  )
})


testthat::test_that("nevent: reproduces documented example values", {
  # Rd example: piecewise accrual, piecewise exponential survivals, 5% dropout.
  out <- nevent(
    time = c(9, 24),
    allocationRatioPlanned = 1,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 12,
    maxFollowupTime = 30
  )

  testthat::expect_equal(
    out,
    matrix(c(13.109896137425,  49.619690296716,
             13.4367060727502, 60.8130481814674),
           nrow = 2, ncol = 2),
    tolerance = 1e-8
  )
})


testthat::test_that("invalid time grids and negative rates raise errors", {
  # piecewiseSurvivalTime must start with 0
  testthat::expect_error(
    dtpwexp(q = 1, piecewiseSurvivalTime = c(1, 6),
            lambda = c(0.05, 0.03)),
    "piecewiseSurvivalTime must start with 0"
  )
  testthat::expect_error(
    ptpwexp(q = 1, piecewiseSurvivalTime = c(2, 6),
            lambda = c(0.05, 0.03)),
    "piecewiseSurvivalTime must start with 0"
  )
  testthat::expect_error(
    qtpwexp(p = 0.5, piecewiseSurvivalTime = c(3, 9),
            lambda = c(0.05, 0.03)),
    "piecewiseSurvivalTime must start with 0"
  )

  # piecewiseSurvivalTime should be increasing (no ties, no decreasing)
  testthat::expect_error(
    dtpwexp(q = 1, piecewiseSurvivalTime = c(0, 5, 3),
            lambda = c(0.05, 0.03, 0.02)),
    "piecewiseSurvivalTime should be increasing"
  )
  testthat::expect_error(
    ptpwexp(q = 1, piecewiseSurvivalTime = c(0, 6, 6),
            lambda = c(0.05, 0.03, 0.02)),
    "piecewiseSurvivalTime should be increasing"
  )

  # lambda length must equal piecewiseSurvivalTime length
  testthat::expect_error(
    dtpwexp(q = 1, piecewiseSurvivalTime = c(0, 6),
            lambda = c(0.05, 0.03, 0.02)),
    "lambda and piecewiseSurvivalTime must have the same length"
  )
  testthat::expect_error(
    rtpwexp(n = 5, piecewiseSurvivalTime = c(0, 6),
            lambda = c(0.05)),
    "lambda and piecewiseSurvivalTime must have the same length"
  )

  # lambda must be nonnegative
  testthat::expect_error(
    dtpwexp(q = 1, piecewiseSurvivalTime = c(0, 6),
            lambda = c(-0.05, 0.03)),
    "lambda must be nonnegative"
  )
  testthat::expect_error(
    ptpwexp(q = 1, piecewiseSurvivalTime = c(0, 6),
            lambda = c(0.05, -0.03)),
    "lambda must be nonnegative"
  )
  testthat::expect_error(
    qtpwexp(p = 0.5, piecewiseSurvivalTime = c(0, 6),
            lambda = c(0.05, -0.03)),
    "lambda must be nonnegative"
  )
})
