library(survival)

testthat::test_that("kmest: estimate and standard error", {
  df1 <- kmest(aml, stratum="x", time="time", event="status",
               conftype="none")

  df2 <- summary(survfit(Surv(time, status) ~ x, data=aml,
                         conf.type="none"))

  testthat::expect_equal(df1$survival, df2$surv)
  testthat::expect_equal(df1$stderr, df2$std.err)
})


testthat::test_that("kmest: confidence interval", {
  for (conftype in c("plain", "log", "log-log", "arcsin")) {
    df1 <- kmest(aml, stratum="x", time="time", event="status",
                 conftype=conftype)

    df2 <- summary(survfit(Surv(time, status) ~ x, data=aml,
                           conf.type=conftype))

    testthat::expect_equal(df1$lower, df2$lower)
    testthat::expect_equal(df1$upper, df2$upper)
  }
})


