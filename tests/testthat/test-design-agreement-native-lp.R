testthat::test_that("getDesignAgreement uses the native LP solver", {
  design <- getDesignAgreement(
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

  testthat::expect_equal(design$n, 82)
  testthat::expect_true(all(design$piH0 >= -1e-10))
  testthat::expect_equal(rowSums(design$piH0), c(0.1, 0.2, 0.3, 0.4), tolerance = 1e-10)
  testthat::expect_equal(colSums(design$piH0), c(0.15, 0.2, 0.24, 0.41), tolerance = 1e-10)
})
