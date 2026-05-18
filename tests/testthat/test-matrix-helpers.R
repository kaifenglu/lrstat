testthat::test_that("matrix decomposition reconstructs inputs within tolerance", {
  A <- matrix(c(1, 0, 0, 0,
                1, 2, 0, 0,
                0, 1, 3, 0,
                0, 0, 1, 4), 4, 4)

  sv <- svdcpp(A)

  testthat::expect_equal(
    sv$d,
    c(4.26000667399714, 3.10734855047139, 2.11178458483242, 0.858541685297334),
    tolerance = 1e-6
  )
  recon <- sv$U %*% diag(sv$d) %*% t(sv$V)
  testthat::expect_equal(unname(recon), unname(A), tolerance = 1e-12)
})


testthat::test_that("helper outputs agree with base R equivalents on fixtures", {
  # float_to_fraction fixtures
  testthat::expect_equal(float_to_fraction(5 / 3), c(5, 3))
  testthat::expect_equal(float_to_fraction(0.125), c(1, 8))
  testthat::expect_equal(float_to_fraction(pi, tol = 1e-5), c(355, 113))

  # fquantile agrees with qweibull when S is Weibull survival.
  probs <- c(0.25, 0.5, 0.75)
  fq <- fquantile(pweibull, probs = probs, shape = 1.37, scale = 1 / 0.818, lower.tail = FALSE)
  qref <- qweibull(probs, shape = 1.37, scale = 1 / 0.818)
  testthat::expect_equal(fq, c(0.492372252261943, 0.935536878350996, 1.55164302597529), tolerance = 1e-6)
  testthat::expect_equal(fq, qref, tolerance = 5e-6)

  # Multiplicity helpers: numeric baselines from documented examples.
  pvalues <- matrix(c(0.01, 0.005, 0.015, 0.022,
                      0.02, 0.015, 0.010, 0.023),
                    nrow = 2, ncol = 4, byrow = TRUE)

  testthat::expect_equal(
    ftrunc(pvalues, "hochberg"),
    matrix(c(0.022, 0.020, 0.022, 0.022,
             0.023, 0.023, 0.023, 0.023), nrow = 2, byrow = TRUE),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    ftrunc(pvalues[1, ], "holm", gamma = 0.8),
    c(0.0315789473684211, 0.02, 0.0333333333333333, 0.0333333333333333),
    tolerance = 1e-12
  )

  p <- c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
  family <- matrix(c(1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 0, 0,
                     0, 0, 0, 0, 0, 0, 1, 1), nrow = 4, byrow = TRUE)
  serial <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0,
                     1, 0, 0, 0, 0, 0, 0, 0,
                     0, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0, 0, 0, 0,
                     0, 0, 0, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 0, 0, 0,
                     0, 0, 0, 0, 0, 1, 0, 0), nrow = 8, byrow = TRUE)
  parallel <- matrix(0, 8, 8)
  gamma <- c(0.6, 0.6, 0.6, 1)

  testthat::expect_equal(
    fstdmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = FALSE),
    c(0.02425, 0.0136, 0.033875, 0.02425, 0.04625, 0.02425, 0.0814, 0.033875),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    fmodmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = TRUE),
    c(0.02425, 0.0136, 0.033, 0.02425, 0.037, 0.02425, 0.0814, 0.033),
    tolerance = 1e-10
  )

  w <- c(0.5, 0.5, 0, 0)
  g <- matrix(c(0, 0, 1, 0,
                0, 0, 0, 1,
                0, 1, 0, 0,
                1, 0, 0, 0), nrow = 4, byrow = TRUE)
  testthat::expect_equal(
    fwgtmat(w, g),
    matrix(c(0.5, 0.5, 0.0, 0.0,
             0.5, 0.5, 0.0, 0.0,
             0.5, 0.5, 0.0, 0.0,
             0.5, 0.5, 0.0, 0.0,
             0.5, 0.0, 0.0, 0.5,
             1.0, 0.0, 0.0, 0.0,
             0.5, 0.0, 0.0, 0.5,
             1.0, 0.0, 0.0, 0.0,
             0.0, 0.5, 0.5, 0.0,
             0.0, 0.5, 0.5, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 0.5, 0.5,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0), nrow = 15, byrow = TRUE),
    tolerance = 1e-12
  )
})


testthat::test_that("invalid dimensions and NA input handling checks", {
  testthat::expect_error(svdcpp(c(1, 2, 3)))

  # Current implementation returns denominator 2 for NA input.
  testthat::expect_equal(float_to_fraction(NA_real_), c(NA_real_, 2))

  testthat::expect_error(
    fquantile(pweibull, probs = 1.1, shape = 1.37, scale = 1 / 0.818, lower.tail = FALSE),
    "probs must be less than or equal to"
  )

  testthat::expect_error(ftrunc(c(0.1, 0.2), test = "bad"), "test must be Holm, Hochberg, or Hommel")
  testthat::expect_error(ftrunc(c(0.1, 0.2), test = "holm", gamma = 1.2), "gamma must be within \\[0, 1\\]")

  testthat::expect_error(
    fmodmix(
      p = c(0.1, 0.2),
      family = matrix(c(1, 0, 0, 1), nrow = 2),
      serial = matrix(0, 2, 2),
      parallel = matrix(0, 2, 2),
      gamma = c(0.5),
      test = "bad"
    ),
    "test must be either Holm, Hochberg, or Hommel"
  )
})
