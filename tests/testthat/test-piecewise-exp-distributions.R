testthat::test_that("d/p/q consistency and inverse checks across cut boundaries", {
  piecewise_survival_time <- c(0, 6, 9, 15)
  lambda <- c(0.025, 0.04, 0.015, 0.007)
  lower_bound <- 4
  q <- c(4, 5.999999, 6, 6.000001, 8, 8.999999, 9, 14.999999, 15, 18)

  density <- dtpwexp(q = q, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound)
  density_log <- dtpwexp(q = q, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound, log.d = TRUE)
  cdf <- ptpwexp(q = q, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound)
  surv <- ptpwexp(q = q, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound, lower.tail = FALSE)
  q_from_cdf <- qtpwexp(p = cdf, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound)
  q_from_surv <- qtpwexp(p = surv, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound, lower.tail = FALSE)
  cdf_log <- ptpwexp(q = q, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda, lowerBound = lower_bound, log.p = TRUE)

  testthat::expect_true(all(is.finite(density)))
  testthat::expect_true(all(density >= 0))
  testthat::expect_equal(density[1], 0)
  testthat::expect_equal(exp(density_log), density, tolerance = 1e-12)
  testthat::expect_equal(cdf + surv, rep(1, length(q)), tolerance = 1e-12)
  testthat::expect_equal(exp(cdf_log), cdf, tolerance = 1e-12)
  testthat::expect_equal(q_from_cdf, q, tolerance = 1e-8)
  testthat::expect_equal(q_from_surv, q, tolerance = 1e-8)
  testthat::expect_true(all(diff(cdf) >= -1e-12))
  testthat::expect_true(all(diff(surv) <= 1e-12))
  testthat::expect_equal(cdf[1], 0)
  testthat::expect_equal(surv[1], 1)
})


testthat::test_that("random generation moments match theoretical moments", {
  piecewise_survival_time <- c(0, 6, 9, 15)
  lambda <- c(0.025, 0.04, 0.015, 0.007)
  theoretical <- mtpwexp(piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda)

  set.seed(20260518)
  x <- rtpwexp(50000, piecewiseSurvivalTime = piecewise_survival_time, lambda = lambda)

  testthat::expect_true(all(is.finite(x)))
  testthat::expect_true(all(x >= 0))
  testthat::expect_equal(mean(x), theoretical$mean, tolerance = 1.5)
  testthat::expect_equal(stats::var(x), theoretical$variance, tolerance = 1000)
})


testthat::test_that("pwexpcuts and pwexploglik behave locally around the true parameters", {
  piecewise_fit <- pwexpcuts(
    ptpwexp,
    piecewiseSurvivalTime = c(0, 3.4, 5.5),
    lambda = c(0.0168, 0.0833, 0.0431),
    lowerBound = 0,
    lower.tail = FALSE
  )
  weibull_fit <- pwexpcuts(
    pweibull,
    shape = 1.37,
    scale = 1 / 0.818,
    lower.tail = FALSE
  )

  testthat::expect_equal(round(piecewise_fit$piecewiseSurvivalTime, 6), c(0, 3.399971, 5.500012))
  testthat::expect_equal(round(piecewise_fit$lambda, 6), c(0.0168, 0.083299, 0.0431))
  testthat::expect_equal(round(piecewise_fit$loglik, 6), c(-4.120403, -4.096664))

  testthat::expect_equal(round(weibull_fit$piecewiseSurvivalTime, 6), c(0, 0.526637, 1.281462, 2.076264, 2.864521, 3.633859))
  testthat::expect_equal(round(weibull_fit$lambda, 6), c(0.584284, 0.975661, 1.23927, 1.43633, 1.593907, 1.762896))
  testthat::expect_true(all(diff(weibull_fit$loglik) > 0))

  synthetic_survival <- function(t) {
    ptpwexp(
      t,
      piecewiseSurvivalTime = c(0, 3.4, 5.5),
      lambda = c(0.0168, 0.0833, 0.0431),
      lower.tail = FALSE
    )
  }

  true_fit <- pwexploglik(c(3.4, 5.5), synthetic_survival)
  left_shift <- pwexploglik(c(3.2, 5.5), synthetic_survival)
  right_shift <- pwexploglik(c(3.6, 5.5), synthetic_survival)
  inner_shift <- pwexploglik(c(3.4, 5.3), synthetic_survival)
  outer_shift <- pwexploglik(c(3.4, 5.7), synthetic_survival)

  testthat::expect_equal(round(true_fit$piecewiseSurvivalTime, 6), c(0, 3.4, 5.5))
  testthat::expect_equal(round(true_fit$lambda, 6), c(0.0168, 0.0833, 0.0431))
  testthat::expect_equal(round(true_fit$loglik, 6), -4.096663)
  testthat::expect_gt(true_fit$loglik, left_shift$loglik)
  testthat::expect_gt(true_fit$loglik, right_shift$loglik)
  testthat::expect_gt(true_fit$loglik, inner_shift$loglik)
  testthat::expect_gt(true_fit$loglik, outer_shift$loglik)
})