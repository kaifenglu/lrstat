library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("phregr: handling ties", {
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

