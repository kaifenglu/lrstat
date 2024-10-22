testthat::test_that("logisregr: freq works", {
  fit1 <- logisregr(ingots, event="NotReady", covariates="Heat*Soak",
                    freq="Freq")

  # expand data into individual subjects
  ingots2 <- ingots[rep(1:nrow(ingots), ingots$Freq),]
  fit2 <- logisregr(ingots2, event="NotReady", covariates="Heat*Soak")

  testthat::expect_equal(fit1$sumstat, fit2$sumstat)
  testthat::expect_equal(fit1$parest, fit2$parest)
})


testthat::test_that("logisregr: firth with plci", {
  fit1 <- logisregr(sexagg, event="case",
                    covariates=c("age", "oc", "vic", "vicl", "vis", "dia"),
                    freq="COUNT", firth=TRUE, plci=TRUE)

  # the following values are almost identical to those from SAS PROC LOGISTIC
  # with numerical differences due to different convergence criteria used
  coef = c(0.1203, -1.1060, -0.0688, 2.2689, -2.1114, -0.7883, 3.0959)
  ci.lower = c(-0.8186, -1.9738, -0.9414, 1.2730, -3.2609, -1.6081, 0.7746)
  ci.upper = c(1.0732, -0.3074, 0.7892, 3.4354, -1.1177, 0.0152, 8.0303)

  testthat::expect_equal(round(fit1$parest$beta, 4), coef)
  testthat::expect_equal(round(fit1$parest$lower, 4), ci.lower)
  testthat::expect_equal(round(fit1$parest$upper, 4), ci.upper)
})


testthat::test_that("logisregr: firth with flic", {
  fit1 <- logisregr(sexagg, event="case",
                    covariates=c("age", "oc", "vic", "vicl", "vis", "dia"),
                    freq="COUNT", firth=TRUE, flic=TRUE)

  # the following value is from the flic function in the logistf package
  # of note, the logistf function with the flic option yields a different
  # value, and the variance of the intercept is off in the logistf package
  intercept = 0.1315

  testthat::expect_equal(round(fit1$parest$beta[1], 4), 0.1315)
})


testthat::test_that("logisregr: robust variance", {
  options(contrasts = c("contr.SAS", "contr.poly"))

  fit1 <- logisregr(six, event="wheeze",
                    covariates=c("city", "age", "smoke"),
                    id="case", robust=TRUE)

  # the following values can be obtained from SAS PROC GENMOD statements:
  #   proc genmod data=six;
  #     class case city;
  #     model wheeze(event='1') = city age smoke / dist=bin;
  #     repeated subject=case / type=ind;
  #   run;
  # here the standard error is based on the robust variance estimate
  coef = c(1.2597, 0.1391, -0.2003, -0.1284)
  stderr = c(3.0645, 0.6859, 0.2820, 0.3926)

  testthat::expect_equal(round(fit1$parest$beta, 4), coef)
  testthat::expect_equal(round(fit1$parest$sebeta, 4), stderr)
})
