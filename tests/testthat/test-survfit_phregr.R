library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that("survfit_phregr: right-censored data", {
  pbc <- pbc %>% mutate(event = 1*(status == 2))

  # fit a Cox model to the original data set
  fit1 <- phregr(pbc, time="time", event="event",
                 covariates=c("age", "edema", "log(bili)",
                              "log(protime)", "log(albumin)"))

  fit2 <- coxph(Surv(time, event) ~ age + edema + log(bili) +
                  log(protime) + log(albumin), data=pbc)

  # create a data set corresponding to the hypothetical subject
  temp <- data.frame(age=53, edema=0, bili=2, protime=12, albumin=2)

  # obtain the expected survival curve
  surv1 <- survfit_phregr(fit1, newdata=temp)
  surv2 <- survfit(fit2, newdata=temp, conf.type="log-log")

  # extract common variables
  surv1b <- surv1 %>%
    select(time, nrisk, nevent, ncensor, cumhaz,
           surv, sesurv, lower, upper) %>%
    slice(-1)

  surv2b <- data.frame(time = surv2$time, nrisk = surv2$n.risk,
                       nevent = surv2$n.event, ncensor = surv2$n.censor,
                       cumhaz = surv2$cumhaz, surv = surv2$surv,
                       sesurv = surv2$surv*surv2$std.chaz,
                       lower = surv2$lower, upper = surv2$upper)

  testthat::expect_equal(surv1b, surv2b)
})


testthat::test_that("survfit_phregr: stratified analysis", {
  pbc <- pbc %>% mutate(event = 1*(status == 2))

  fit1 <- phregr(pbc, time="time", event="event",
                 covariates=c("age", "log(bili)"), stratum="edema")

  fit2 <- coxph(Surv(time, event) ~ age + log(bili) + strata(edema),
                data=pbc)

  # create a data set corresponding to the hypothetical subject
  # of note, survfit_phregr requests explict stratum info in newdata
  # while survfit.coxph replicates the covariate values for each stratum
  temp1 <- data.frame(edema=c(0, 0.5, 1)) %>%
    cross_join(data.frame(age=c(53, 60), bili = c(2,3)))
  temp2 <- data.frame(age=c(53, 60), bili = c(2,3))

  # obtain the expected survival curve
  surv1 <- survfit_phregr(fit1, newdata=temp1)
  surv2 <- survfit(fit2, newdata=temp2, conf.type="log-log")

  # of note, surv1 is one data frame ordered by strata and covariates,
  # in contrast, surv2 has the strata in rows and covariates in columns
  # need to rearrange surv2 to have the same layout as surv1
  surv1b <- surv1 %>%
    mutate(bili = exp(log.bili.)) %>%
    select(time, nrisk, nevent, cumhaz, surv,
           sesurv, lower, upper, edema, age, bili) %>%
    group_by(edema, age, bili) %>%
    slice(-1) %>%
    ungroup()

  strata <- rep(as.numeric(substring(names(surv2$strata), 7)), surv2$strata)

  surv2b <- bind_rows(
    tibble(time = surv2$time, nrisk = surv2$n.risk,
           nevent = surv2$n.event, cumhaz = surv2$cumhaz[,1],
           surv = surv2$surv[,1],
           sesurv = surv2$surv[,1]*surv2$std.chaz[,1],
           lower = surv2$lower[,1], upper = surv2$upper[,1],
           edema = strata) %>% bind_cols(surv2$newdata[1,]),
    tibble(time = surv2$time, nrisk = surv2$n.risk,
           nevent = surv2$n.event, cumhaz = surv2$cumhaz[,2],
           surv = surv2$surv[,2],
           sesurv = surv2$surv[,2]*surv2$std.chaz[,2],
           lower = surv2$lower[,2], upper = surv2$upper[,2],
           edema = strata) %>%
      bind_cols(surv2$newdata[2,])) %>%
    arrange(edema)

  testthat::expect_equal(surv1b, surv2b)
})


testthat::test_that("survfit_phregr: time-dependent covariates", {
  # create counting process data
  pbcseq <- pbcseq %>%
    group_by(id) %>%
    arrange(id, day) %>%
    mutate(tstart = day,
           tstop = ifelse(row_number() != n(), lead(day), futime),
           event = ifelse(row_number() != n(), 0, 1*(status == 2)))

  # fit a Cox model to the original data, note the use of robust variance
  fit1 <- phregr(pbcseq, time="tstart", time2="tstop", event="event",
                 covariates=c("age", "edema", "log(bili)",
                              "log(protime)", "log(albumin)"),
                 id="id", robust=TRUE)

  fit2 <- coxph(Surv(tstart, tstop, event) ~ age + edema + log(bili) +
                  log(protime) + log(albumin), data=pbcseq,
                cluster=id, robust=TRUE)

  # create a data set corresponding to the hypothetical subject with
  # time-dependent covariates
  temp <- data.frame(id      = c(   0,    0,    0,    0,    0),
                     tstart  = c(   0,  365,  730, 1095, 1460),
                     tstop   = c( 365,  730, 1095, 1460, 3000),
                     event   = c(   0,    0,    0,    0,    0),
                     age     = c(  53,   53,   53,   53,   53),
                     bili    = c(   1,    2,    3,    5,    7),
                     edema   = c(   1,    1,    1,    1,    1),
                     albumin = c( 3.5,  3.5,  3.5,  3.5,  3.5),
                     protime = c(  11,   11,   11,   11,   11))

  surv1 <- survfit_phregr(fit1, newdata=temp)
  surv2 <- survfit(fit2, newdata=temp, id=id, conf.type="log-log")

  surv1b <- surv1 %>%
    select(time, nrisk, nevent, ncensor, cumhaz,
           surv, sesurv, lower, upper) %>%
    slice(-1)

  surv2b <- data.frame(time = surv2$time, nrisk = surv2$n.risk,
                       nevent = surv2$n.event, ncensor = surv2$n.censor,
                       cumhaz = surv2$cumhaz, surv = surv2$surv,
                       sesurv = surv2$surv*surv2$std.chaz,
                       lower = surv2$lower, upper = surv2$upper)

  # the point estimate would match but the standard error and confidence
  # interval do not match due to an apparent bug in the survival package
  testthat::expect_equal(surv1b %>% select(time, surv),
                         surv2b %>% select(time, surv))
})

