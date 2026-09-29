testthat::test_that("lrsim_pkTrtSel supports all Carreras selection rules", {
  for (rule in c("exposure", "os", "perfect")) {
    simulation <- lrsim_pkTrtSel(
      selectionRule = rule,
      maxNumberOfIterations = 20L,
      seed = 314159L,
      nthreads = 1L
    )

    testthat::expect_s3_class(simulation, "lrsim_pkTrtSel")
    testthat::expect_equal(simulation$numberOfIterations, 20)
    testthat::expect_equal(simulation$stage1Patients,
                 simulation$interimPatients)
    testthat::expect_equal(nrow(simulation$sumdata), 20)
    testthat::expect_setequal(
      simulation$methods,
      c("followupwise", "patientwise", "dunnett", "stage2only", "threearm")
    )
    probabilities <- vapply(
      simulation$byMethod,
      function(method) method$rejectionProbability,
      numeric(1)
    )
    testthat::expect_true(all(probabilities >= 0 & probabilities <= 1))
  }
})

testthat::test_that("lrsim_pkTrtSel summaries reproduce reported results", {
  simulation <- lrsim_pkTrtSel(
    medianSurvival = c(9, 7.5, 6),
    meanExposure = c(300, 185),
    selectionRule = "exposure",
    maxNumberOfIterations = 40L,
    seed = 271828L,
    nthreads = 1L
  )

  testthat::expect_equal(
    simulation$probabilitySelectHigh,
    mean(simulation$sumdata$selectedRegimen == 1)
  )
  testthat::expect_equal(
    simulation$probabilityCorrectSelection,
    mean(simulation$sumdata$correctSelection)
  )
  for (method in simulation$methods) {
    testthat::expect_equal(
      simulation$byMethod[[method]]$rejectionProbability,
      mean(simulation$sumdata[[paste0("reject.", method)]])
    )
  }
})

testthat::test_that("lrsim_pkTrtSel validates article design inputs", {
  testthat::expect_error(
    lrsim_pkTrtSel(medianSurvival = c(9, 6)),
    "medianSurvival"
  )
  testthat::expect_error(
    lrsim_pkTrtSel(corrSurvivalExposure = 1.1),
    "corrSurvivalExposure"
  )
})

testthat::test_that("Carreras treatment-selection benchmarks are reproduced", {
  exposure <- lrsim_pkTrtSel(
    medianSurvival = c(9, 7.5, 6),
    meanExposure = c(300, 175),
    exposureSD = 120,
    corrSurvivalExposure = 0.2,
    interimPatients = 100L,
    selectionRule = "exposure",
    maxNumberOfIterations = 5000L,
    seed = 314159L,
    nthreads = 1L
  )
  testthat::expect_equal(exposure$probabilitySelectHigh, 0.9, tolerance = 0.03)

  perfect <- lrsim_pkTrtSel(
    medianSurvival = c(6, 6, 6),
    meanExposure = c(300, 200),
    exposureSD = 120,
    corrSurvivalExposure = 0,
    interimPatients = 100L,
    selectionRule = "perfect",
    maxNumberOfIterations = 5000L,
    seed = 271828L,
    nthreads = 1L
  )
  testthat::expect_gt(
    perfect$byMethod$followupwise$rejectionProbability,
    0.025
  )
})