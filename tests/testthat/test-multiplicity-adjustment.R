testthat::test_that("error spending is non-decreasing and ends at the target error", {
  # Single look: sfOF at t=0.5 gives a small fraction of the total error.
  es_single <- errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
  testthat::expect_true(es_single > 0 && es_single < 0.025)
  testthat::expect_equal(es_single, 0.001525323, tolerance = 1e-6)

  # Three-look sfOF: strictly non-decreasing, terminates exactly at alpha.
  es_of <- errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfOF")
  testthat::expect_equal(length(es_of), 3L)
  testthat::expect_true(all(diff(es_of) > 0))
  testthat::expect_equal(es_of[3], 0.025)
  testthat::expect_equal(es_of, c(0.001525323, 0.009649325, 0.025), tolerance = 1e-6)

  # sfHSD with gamma=-4 also ends at alpha and is non-decreasing.
  es_hsd <- errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
  testthat::expect_true(all(diff(es_hsd) > 0))
  testthat::expect_equal(es_hsd[3], 0.025)
  testthat::expect_equal(es_hsd, c(0.002980073, 0.008902144, 0.025), tolerance = 1e-6)

  # Ten equally-spaced information fractions: differences are all positive.
  t_seq <- seq(0.1, 1, by = 0.1)
  es_seq <- errorSpent(t_seq, error = 0.025, sf = "sfOF")
  testthat::expect_true(all(diff(es_seq) > 0))
  testthat::expect_equal(es_seq[10], 0.025)

  # Pocock spending: gentler at the start, also ends at alpha.
  es_p <- errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfP")
  testthat::expect_true(all(diff(es_p) > 0))
  testthat::expect_equal(es_p[3], 0.025)
  testthat::expect_equal(es_p, c(0.01550286, 0.02069972, 0.025), tolerance = 1e-6)
})


testthat::test_that("Bonferroni, Dunnett, and Simes benchmarks match known references", {
  # Shared setup: two-look pvalues, weight/transition graph.
  pvalues <- matrix(
    c(0.01, 0.005, 0.015, 0.022,
      0.02, 0.015, 0.010, 0.023),
    nrow = 2, ncol = 4, byrow = TRUE
  )
  w <- c(0.5, 0.5, 0, 0)
  G <- matrix(
    c(0, 0, 1, 0,
      0, 0, 0, 1,
      0, 1, 0, 0,
      1, 0, 0, 0),
    nrow = 4, ncol = 4, byrow = TRUE
  )
  wgtmat <- fwgtmat(w, G)
  family <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 2, ncol = 4, byrow = TRUE)
  corr   <- matrix(
    c(1, 0.5, NA, NA,
      0.5, 1, NA, NA,
      NA, NA, 1, 0.5,
      NA, NA, 0.5, 1),
    nrow = 4, byrow = TRUE
  )

  # Bonferroni (fadjpbon).
  bon <- fadjpbon(wgtmat, pvalues)$padj
  testthat::expect_equal(nrow(bon), 2L)
  testthat::expect_equal(ncol(bon), 4L)
  testthat::expect_true(all(bon >= 0 & bon <= 1))
  testthat::expect_equal(
    bon,
    matrix(c(0.02, 0.01, 0.03, 0.03,
             0.04, 0.03, 0.04, 0.04),
           nrow = 2, ncol = 4, byrow = TRUE),
    tolerance = 1e-10
  )
  # Adjusted p-values must be >= the corresponding raw p-values.
  testthat::expect_true(all(bon >= pvalues))

  # Dunnett (fadjpdun).
  dun <- fadjpdun(wgtmat, pvalues, family, corr)$padj
  testthat::expect_equal(dim(dun), c(2L, 4L))
  testthat::expect_true(all(dun >= 0 & dun <= 1))
  testthat::expect_equal(
    dun,
    matrix(c(0.02,       0.01,       0.02772937, 0.02772937,
             0.04,       0.02772937, 0.04,       0.04),
           nrow = 2, ncol = 4, byrow = TRUE),
    tolerance = 1e-6
  )
  testthat::expect_true(all(dun >= pvalues))

  # Simes (fadjpsim).
  sim <- fadjpsim(wgtmat, pvalues, family)$padj
  testthat::expect_equal(dim(sim), c(2L, 4L))
  testthat::expect_true(all(sim >= 0 & sim <= 1))
  testthat::expect_equal(
    sim,
    matrix(c(0.02, 0.01, 0.022, 0.022,
             0.04, 0.02, 0.04,  0.04),
           nrow = 2, ncol = 4, byrow = TRUE),
    tolerance = 1e-10
  )
  testthat::expect_true(all(sim >= pvalues))

  # Dunnett is at least as powerful as Bonferroni (adjusted p-values are <=).
  testthat::expect_true(all(dun <= bon + 1e-10))
})


testthat::test_that("sequential procedures produce valid outcomes", {
  # fseqbon: group sequential Bonferroni — Maurer & Bretz (2013) case study.
  fseq_out <- fseqbon(
    w = c(0.5, 0.5, 0, 0),
    G = matrix(
      c(0, 0.5, 0.5, 0,
        0.5, 0, 0, 0.5,
        0, 1, 0, 0,
        1, 0, 0, 0),
      nrow = 4, ncol = 4, byrow = TRUE
    ),
    alpha = 0.025,
    kMax = 3,
    typeAlphaSpending = rep("sfOF", 4),
    maxInformation = rep(1, 4),
    k1 = 2,
    p = matrix(
      c(0.0062, 0.017, 0.009, 0.13,
        0.0002, 0.0035, 0.002, 0.06),
      nrow = 2, ncol = 4, byrow = TRUE
    ),
    information = matrix(
      c(rep(1/3, 4), rep(2/3, 4)),
      nrow = 2, ncol = 4, byrow = TRUE
    ),
    nthreads = 1
  )

  # H1, H2, H3 are rejected at look 2; H4 is never rejected.
  testthat::expect_equal(fseq_out, c(2, 2, 2, 0))
  # All values must be in {0, 1, ..., kMax}.
  testthat::expect_true(all(fseq_out >= 0 & fseq_out <= 3))

  # fstp2seq: stepwise gatekeeping — Hochberg with retesting.
  p_raw <- c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
  gamma <- c(0.6, 0.6, 0.6, 1)
  stp2 <- fstp2seq(p_raw, gamma, test = "hochberg", retest = 1)

  testthat::expect_equal(length(stp2), 8L)
  testthat::expect_true(all(stp2 >= 0 & stp2 <= 1))
  testthat::expect_equal(
    stp2,
    c(0.02425, 0.01360, 0.03300, 0.02425, 0.03700, 0.02425, 0.08140, 0.03300),
    tolerance = 1e-5
  )
  # All adjusted p-values must be >= the corresponding raw p-values.
  testthat::expect_true(all(stp2 >= p_raw - 1e-10))

  # repeatedPValue: single-sequence example — p-values decrease over looks.
  rp1 <- repeatedPValue(
    kMax = 3,
    typeAlphaSpending = "sfOF",
    maxInformation = 800,
    p = c(0.2, 0.15, 0.1),
    information = c(529, 700, 800),
    spendingTime = c(0.6271186, 0.8305085, 1),
    nthreads = 1
  )
  testthat::expect_equal(length(rp1), 3L)
  testthat::expect_true(all(rp1 >= 0 & rp1 <= 1))
  testthat::expect_equal(
    rp1,
    c(0.3101673, 0.2258992, 0.1232186),
    tolerance = 1e-6
  )
  # Repeated p-values are strictly decreasing when raw p-values decrease
  # and information accumulates toward maxInformation.
  testthat::expect_true(all(diff(rp1) < 0))

  # repeatedPValue: multi-hypothesis matrix input.
  rp2 <- repeatedPValue(
    kMax = 3,
    typeAlphaSpending = "sfOF",
    p = matrix(
      c(0.0062, 0.017, 0.009, 0.13,
        0.0002, 0.0035, 0.002, 0.06),
      nrow = 4, ncol = 2
    ),
    information = c(1/3, 2/3),
    nthreads = 1
  )
  testthat::expect_equal(dim(rp2), c(4L, 2L))
  testthat::expect_true(all(rp2 >= 0 & rp2 <= 1))
  testthat::expect_equal(
    rp2,
    matrix(
      c(0.1140577, 0.002393359,
        0.1682137, 0.017160592,
        0.1315366, 0.011648306,
        0.3820275, 0.128460197),
      nrow = 4, ncol = 2, byrow = TRUE
    ),
    tolerance = 1e-5
  )
})


testthat::test_that("fseqbon lookback handling depends on visit schedule and alpha recycling", {
  incidence_matrix <- matrix(
    c(1, 0, 0, 0, 0,
      0, 1, 1, 0, 0,
      0, 0, 1, 1, 1),
    nrow = 5,
    ncol = 3
  )
  w <- c(0.001, 0.024, 0) / 0.025
  G <- matrix(
    c(0, 1, 0,
      0, 0, 1,
      1, 0, 0),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )

  max_information <- c(600, 120, 100)

  information <- matrix(
    c(600, 24, 8,
      600, 54, 24,
      600, 120, 48,
      600, 120, 72,
      600, 120, 100),
    nrow = 5,
    byrow = TRUE
  )

  # Scenario 1: PFS and OS p-values decrease over looks. OS is first
  # rejected at the final look, and ORR can be rejected only via lookback
  # after late alpha recycling from OS.
  pvals_lookback_sensitive <- matrix(
    c(0.006, 0.9, 0.9,
      0.006, 0.006, 0.9,
      0.006, 0.003, 0.03,
      0.006, 0.003, 0.02,
      0.006, 0.003, 0.001),
    nrow = 5,
    byrow = TRUE
  )

  out_true <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfP", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_lookback_sensitive,
    information = information,
    lookback = TRUE,
    nthreads = 1
  )
  out_false <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfP", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_lookback_sensitive,
    information = information,
    lookback = FALSE,
    nthreads = 1
  )

  testthat::expect_equal(out_true, c(1, 2, 5))
  testthat::expect_equal(out_false, c(0, 2, 5))

  # Scenario 2: ORR is already significant at look 1, so lookback toggle
  # should not change the rejection pattern.
  pvals_not_sensitive <- pvals_lookback_sensitive
  pvals_not_sensitive[1, 1] <- 0.0008

  out_true_2 <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfP", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_not_sensitive,
    information = information,
    lookback = TRUE,
    nthreads = 1
  )
  out_false_2 <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfP", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_not_sensitive,
    information = information,
    lookback = FALSE,
    nthreads = 1
  )

  testthat::expect_equal(out_true_2, c(1, 2, 5))
  testthat::expect_equal(out_false_2, c(1, 2, 5))

  # Scenario 3: sfOF spending with decreasing p-values over time for each
  # tested endpoint, while retaining lookback sensitivity.
  pvals_of_lookback_sensitive <- matrix(
    c(0.0035, 0.9, 0.9,
      0.0035, 0.0045, 0.9,
      0.0035, 0.0020, 0.0030,
      0.0035, 0.0020, 0.0015,
      0.0035, 0.0020, 0.0005),
    nrow = 5,
    byrow = TRUE
  )

  out_true_3 <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfOF", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_of_lookback_sensitive,
    information = information,
    lookback = TRUE,
    nthreads = 1
  )
  out_false_3 <- fseqbon(
    w = w,
    G = G,
    alpha = 0.025,
    kMax = 5,
    typeAlphaSpending = rep("sfOF", 3),
    maxInformation = max_information,
    incidenceMatrix = incidence_matrix,
    k1 = 5,
    p = pvals_of_lookback_sensitive,
    information = information,
    lookback = FALSE,
    nthreads = 1
  )

  testthat::expect_equal(out_true_3, c(1, 3, 3))
  testthat::expect_equal(out_false_3, c(0, 3, 3))
})
