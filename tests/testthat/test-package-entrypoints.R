testthat::test_that("package-level wrappers return expected object types", {
  lr1 <- lrstat(
    time = c(22, 40),
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE,
    predictTarget = 1
  )

  testthat::expect_s3_class(lr1, "data.frame")
  testthat::expect_equal(
    names(lr1),
    c("time", "subjects", "nevents", "nevents1", "nevents2",
      "ndropouts", "ndropouts1", "ndropouts2", "nfmax", "nfmax1", "nfmax2")
  )
  testthat::expect_equal(nrow(lr1), 2L)
  testthat::expect_equal(lr1$subjects, c(468, 468))
  testthat::expect_equal(lr1$nevents, c(154.403999115898, 307.50780621924), tolerance = 1e-5)
})


testthat::test_that("shiny app launcher resolves directory and returns shiny app object or launch handle", {
  testthat::skip_if_not_installed("shiny")
  testthat::skip_if_not_installed("shinyjs")
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tidyr")
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("plotly")

  app_dir <- system.file("shinyApp", package = "lrstat")
  testthat::expect_true(nzchar(app_dir))
  testthat::expect_true(dir.exists(app_dir))
  testthat::expect_true(file.exists(file.path(app_dir, "app.R")))

  app <- runShinyApp_lrstat()
  testthat::expect_s3_class(app, "shiny.appobj")
  testthat::expect_true(is.list(app))
  testthat::expect_true(all(c("httpHandler", "serverFuncSource", "options") %in% names(app)))
})


testthat::test_that("updateGraph validates inputs and returns plot object contracts", {
  ug <- updateGraph(
    w = c(0.5, 0.5, 0, 0),
    G = matrix(c(0, 0.5, 0.5, 0,
                 0.5, 0, 0, 0.5,
                 0, 1, 0, 0,
                 1, 0, 0, 0),
               nrow = 4, ncol = 4, byrow = TRUE),
    I = c(1, 2, 3, 4),
    j = 1
  )

  testthat::expect_type(ug, "list")
  testthat::expect_named(ug, c("w", "G", "I"))
  testthat::expect_equal(ug$w, c(0, 0.75, 0.25, 0), tolerance = 1e-12)
  testthat::expect_equal(ug$I, c(2, 3, 4))
  testthat::expect_equal(
    ug$G,
    matrix(c(0, 0, 0, 0,
             0, 0, 1/3, 2/3,
             0, 1, 0, 0,
             0, 0.5, 0.5, 0), nrow = 4, byrow = TRUE),
    tolerance = 1e-6
  )

  lr2 <- lrstat(
    time = c(22, 40),
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE,
    predictTarget = 2
  )
  testthat::expect_equal(lr2$uscore, c(-6.40409666092687, -24.3518052400386), tolerance = 1e-5)
  testthat::expect_equal(lr2$vscore, c(38.5649730098823, 76.1893888757008), tolerance = 1e-5)
  testthat::expect_equal(lr2$logRankZ, c(-1.03124408865089, -2.78986991074603), tolerance = 1e-5)
  testthat::expect_equal(lr2$hazardRatioH0, c(1, 1))

  lr3 <- lrstat(
    time = c(22, 40),
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE,
    predictTarget = 3
  )
  testthat::expect_equal(lr3$HR, c(0.846828140276223, 0.726621845196838), tolerance = 1e-5)
  testthat::expect_equal(lr3$vlogHR, c(0.0260517943066541, 0.0132126583490171), tolerance = 1e-5)
  testthat::expect_equal(lr3$zlogHR, c(-1.03005961198032, -2.77824545182258), tolerance = 1e-5)
})
