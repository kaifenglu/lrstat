testthat::test_that("probabilities are in [0, 1] and monotone in bounds", {
  p_rect <- pbvnorm(c(-1, -1), c(1, 1), 0.5)
  p_quadrant_ind <- pbvnorm(c(-Inf, -Inf), c(0, 0), 0)
  p_quadrant_pos <- pbvnorm(c(-Inf, -Inf), c(0, 0), 0.7)

  testthat::expect_equal(p_rect, 0.497971777839208, tolerance = 1e-10)
  testthat::expect_equal(p_quadrant_ind, 0.25, tolerance = 1e-12)
  testthat::expect_equal(p_quadrant_pos, 0.373408344446682, tolerance = 1e-10)

  n <- 5
  lower <- rep(-1, n)
  upper <- rep(3, n)
  mean <- rep(0, n)
  sigma_cs <- matrix(0.5, n, n)
  diag(sigma_cs) <- 1

  p1 <- pmvnormr(lower, upper, mean, sigma_cs)
  p2 <- pmvnormr(lower, upper + 0.2, mean, sigma_cs)

  testthat::expect_true(p1 >= 0 && p1 <= 1)
  testthat::expect_true(p2 >= 0 && p2 <= 1)
  testthat::expect_gt(p2, p1)
  testthat::expect_equal(as.numeric(p1), 0.580047666534787, tolerance = 1e-6)
  testthat::expect_equal(attr(p1, "method"), "analytic")
  testthat::expect_equal(attr(p1, "error"), 0)
  testthat::expect_equal(attr(p1, "nsamples"), 1)
})


testthat::test_that("q and p inversion holds within tolerance", {
  n <- 5
  mean <- rep(0, n)

  sigma_cs <- matrix(0.5, n, n)
  diag(sigma_cs) <- 1
  q1 <- qmvnormr(0.5, mean = mean, sigma = sigma_cs)
  p_back1 <- pmvnormr(lower = rep(-Inf, n), upper = rep(q1, n), mean = mean, sigma = sigma_cs)

  testthat::expect_equal(q1, 0.815034671372329, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(p_back1), 0.5, tolerance = 1e-8)

  sigma_g <- matrix(c(1, 0.5, 0.3, 0.2, 0.1,
                      0.5, 1, 0.4, 0.3, 0.2,
                      0.3, 0.4, 1, 0.5, 0.3,
                      0.2, 0.3, 0.5, 1, 0.4,
                      0.1, 0.2, 0.3, 0.4, 1), nrow = n)
  q2 <- qmvnormr(0.7, mean = mean, sigma = sigma_g, seed = 314159, parallel = FALSE)
  p_back2 <- pmvnormr(lower = rep(-Inf, n), upper = rep(q2, n), mean = mean, sigma = sigma_g,
                      seed = 314159, parallel = FALSE)

  testthat::expect_equal(q2, 1.35347786951103, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(p_back2), 0.7, tolerance = 2e-4)
})


testthat::test_that("singular or near-singular covariance handling is validated", {
  sigma_sing <- matrix(c(1, 1, 1, 1), nrow = 2)
  p_sing <- pmvnormr(lower = c(-1, -1), upper = c(1, 1), sigma = sigma_sing,
                     seed = 123, parallel = FALSE)

  testthat::expect_true(is.finite(p_sing))
  testthat::expect_true(p_sing >= 0 && p_sing <= 1)
  testthat::expect_equal(as.numeric(p_sing), 0.341344746068543, tolerance = 1e-6)

  sigma_near <- matrix(c(1, 0.9999, 0.9999, 1), nrow = 2)
  p_near <- pmvnormr(lower = c(-1, -1), upper = c(1, 1), sigma = sigma_near,
                     seed = 123, parallel = FALSE)
  q_near <- qmvnormr(0.5, sigma = sigma_near, seed = 123, parallel = FALSE)

  testthat::expect_true(is.finite(p_near))
  testthat::expect_true(p_near >= 0 && p_near <= 1)
  testthat::expect_equal(as.numeric(p_near), 0.679959105632224, tolerance = 1e-6)
  testthat::expect_equal(attr(p_near, "method"), "analytic")
  testthat::expect_true(is.finite(q_near))
  testthat::expect_equal(q_near, 0.00564187925065574, tolerance = 1e-6)

  testthat::expect_error(pmvnormr(sigma = c(1, 2, 3)), "sigma must be a matrix")
  testthat::expect_error(qmvnormr(0.5, sigma = matrix(1:6, nrow = 2)), "sigma must be square")
})
