testthat::test_that("phregr: handling ties", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  pbc <- pbc %>% mutate(event = 1*(status == 2))

  for (ties in c("breslow", "efron")) {
    fit1 <- phregr(pbc, time="time", event="event",
                   covariates=c("age", "edema", "log(bili)",
                                "log(protime)", "log(albumin)"),
                   ties=ties)

    fit2 <- coxph(Surv(time, event) ~ age + edema + log(bili) +
                    log(protime) + log(albumin), data=pbc, ties=ties)

    dimnames(fit1$vbeta) <- NULL

    testthat::expect_equal(fit1$beta, fit2$coefficients)
    testthat::expect_equal(fit1$vbeta, fit2$var)
  }
})


testthat::test_that("phregr: counting process form", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  heart <- heart %>% mutate(rx = as.numeric(transplant) - 1)

  fit1 <- phregr(heart, time="start", time2="stop", event="event",
                 covariates=c("rx", "age"), id="id", robust=TRUE)

  fit2 <- coxph(Surv(start, stop, event) ~ rx + age, cluster=id,
                data=heart, robust=TRUE)

  dimnames(fit1$vbeta) <- NULL

  testthat::expect_equal(as.numeric(fit1$beta),
                         as.numeric(fit2$coefficients))
  testthat::expect_equal(fit1$vbeta, fit2$var)
  testthat::expect_equal(c(fit1$sumstat$loglik0, fit1$sumstat$loglik1),
                         fit2$loglik)
  testthat::expect_equal(fit1$sumstat$scoretest, fit2$score)
})


testthat::test_that("phregr: firth with plci", {
  library(survival)
  # we include the status variable as a predictor to force an infinite beta
  # in this case, we invoke the firth option to obtain finite beta estimate
  fit1 <- phregr(ovarian, time="futime", event="fustat",
                 covariates=c("rx", "fustat"),
                 firth=TRUE, plci=TRUE)

  # coxph does not have the firth option, and we use SAS PROC PHREG
  #   proc phreg data=ovarian;
  #     model futime*fustat(0) = rx fustat / firth risklimits = pl;
  #   run;
  # to obtain the estimated hazard ratios and 95% profile likelihood CI
  # of note, although not applicable to the ovarian data, which does not
  # have ties, SAS PROC PHREG only has the firth option for the Breslow
  # method for handling ties, while the phregr also has the firth option
  # for the Efron method for handling ties
  beta = c(-0.54197, 4.23615)
  hrlower = c(0.173, 8.771)
  hrupper = c(1.884, 8936.061)

  testthat::expect_equal(round(c(fit1$parest$beta, fit1$parest$lower,
                                 fit1$parest$upper), 3),
                         round(c(beta, log(c(hrlower, hrupper))), 3))
})



testthat::test_that("liferegr: eligible distributions", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  for (dist in c("exponential", "weibull", "lognormal", "loglogistic")) {
    fit1 <- liferegr(ovarian, time="futime", event="fustat",
                     covariates=c("ecog.ps", "rx"), dist=dist)

    fit2 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian,
                    dist=dist)

    testthat::expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
    if (dist != "exponential") {
      testthat::expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
    }
  }
})


testthat::test_that("liferegr: left-censored data", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  tobin <- tobin %>% mutate(time = ifelse(durable > 0, durable, NA))

  for (dist in c("gaussian", "logistic")) {
    fit1 <- liferegr(tobin, time="time", time2="durable",
                     covariates=c("age", "quant"), dist=dist)

    fit2 <- survreg(Surv(durable, durable>0, type='left') ~ age + quant,
                    data=tobin, dist=dist)

    testthat::expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
    testthat::expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
  }
})


testthat::test_that("liferegr: stratification", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  lung1 <- lung %>% mutate(event = 1*(status == 2))

  fit1 <- liferegr(lung1, stratum="sex", time="time", event="event",
                   covariates=c("ph.ecog", "age"))

  fit2 <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex),
                  data=lung)

  testthat::expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
  testthat::expect_equal(exp(fit1$parest$beta[(fit1$nvar+2):(fit1$p)]),
               as.vector(fit2$scale))
})


testthat::test_that("liferegr: robust variance", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
  diabetic <- diabetic %>% mutate(juvenile = 1*(age < 20))

  fit1 <- liferegr(diabetic, time="time", event="status",
                   covariates=c("trt", "juvenile"), id="id",
                   robust=TRUE)
  fit2 <- survreg(Surv(time, status) ~ trt + juvenile, cluster=id,
                  data=diabetic, robust=TRUE)

  testthat::expect_equal(fit1$beta[1:(fit1$nvar+1)], fit2$coefficients)
  testthat::expect_equal(exp(fit1$parest$beta[(fit1$nvar+2)]), fit2$scale)
  testthat::expect_equal(fit1$vbeta, fit2$var)
})



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



# of note, the p-value from lrtest is one-sided pnorm(z), while the p-value
# from survdiff is two-sided 1 - pchisq(z^2, 1)
testthat::test_that("lrtest: two-sided p-value", {
  library(survival)
  df1 <- lrtest(ovarian, treat="rx", time="futime", event="fustat")

  df2 <- survdiff(Surv(futime, fustat) ~ rx, data=ovarian)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})


testthat::test_that("lrtest: stratified log-rank test", {
  library(survival)
  data1 <- subset(rawdata, iterationNumber == 1)

  df1 <- lrtest(data1, stratum="stratum", treat="treatmentGroup",
                time="timeUnderObservation", event="event")

  df2 <- survdiff(Surv(timeUnderObservation, event) ~ treatmentGroup +
                    strata(stratum), data=data1)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})


testthat::test_that("lrtest: weighted log-rank test", {
  library(survival)
  df1 <- lrtest(aml, treat="x", time="time", event="status",
                rho1=0.5)

  df2 <- survdiff(Surv(time, status) ~ x, rho=0.5, data=aml)

  testthat::expect_equal(df1$logRankPValue, df2$pvalue)
})
