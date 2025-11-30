library(survival)

# of note, the p-value from lrtest is one-sided pnorm(z), while the p-value
# from survdiff is two-sided 1 - pchisq(z^2, 1)
testthat::test_that("lrtest: two-sided p-value", {
  df1 <- lrtest(ovarian, treat="rx", time="futime", event="fustat")

  df2 <- survdiff(Surv(futime, fustat) ~ rx, data=ovarian)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})


testthat::test_that("lrtest: stratified log-rank test", {
  data1 <- subset(rawdata, iterationNumber == 1)

  df1 <- lrtest(data1, stratum="stratum", treat="treatmentGroup",
                time="timeUnderObservation", event="event")

  df2 <- survdiff(Surv(timeUnderObservation, event) ~ treatmentGroup +
                    strata(stratum), data=data1)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})


testthat::test_that("lrtest: weighted log-rank test", {
  df1 <- lrtest(aml, treat="x", time="time", event="status",
                rho1=0.5)

  df2 <- survdiff(Surv(time, status) ~ x, rho=0.5, data=aml)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})

