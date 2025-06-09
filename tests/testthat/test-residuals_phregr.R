library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that(
  "residuals_phregr: null model for right-censored data", {
    pbc <- pbc %>% mutate(event = 1*(status == 2))

    for (type in c("martingale", "deviance")) {
      fit1 <- phregr(pbc, time="time", event="event", covariates="")
      rr1 <- residuals_phregr(fit1, type=type)

      fit2 <- coxph(Surv(time, event) ~ 1, data=pbc)
      rr2 <- resid(fit2, type=type)

      names(rr2) <- NULL
      testthat::expect_equal(rr1, rr2)
    }
  })


testthat::test_that(
  "residuals_phregr: right-censored data with covariates", {
    pbc <- pbc %>% mutate(event = 1*(status == 2))

    for (type in c("martingale", "deviance", "score", "schoenfeld",
                   "dfbeta", "dfbetas", "scaledsch")) {
      fit1 <- phregr(pbc, time="time", event="event",
                     covariates=c("age", "edema", "log(bili)",
                                  "log(protime)", "log(albumin)"))

      rr1 <- residuals_phregr(fit1, type=type)

      fit2 <- coxph(Surv(time, event) ~ age + edema + log(bili) +
                      log(protime) + log(albumin), data=pbc)

      rr2 <- resid(fit2, type=type)

      if (type=="martingale" || type=="deviance") {
        names(rr2) <- NULL
      } else {
        rownames(rr2) <- NULL
        if (type=="schoenfeld" || type=="scaledsch") {
          attr(rr1, "time") <- NULL
        } else if (type=="dfbeta" || type=="dfbetas") {
          colnames(rr1) <- NULL
        }
      }

      testthat::expect_equal(rr1, rr2)
    }
  })


testthat::test_that(
  "residuals_phregr: null model for counting process data", {
    pbcseq <- pbcseq %>%
      group_by(id) %>%
      arrange(id, day) %>%
      mutate(tstart = day,
             tstop = ifelse(row_number() != n(), lead(day), futime),
             event = ifelse(row_number() != n(), 0, 1*(status == 2)))

    for (type in c("martingale", "deviance")) {
      fit1 <- phregr(pbcseq, time="tstart", time2="tstop", event="event",
                     covariates="", id="id", robust=TRUE)

      rr1 <- residuals_phregr(fit1, type=type)

      fit2 <- coxph(Surv(tstart, tstop, event) ~ 1, data=pbcseq,
                    cluster=id, robust=TRUE)

      rr2 <- resid(fit2, type=type)

      names(rr2) <- NULL
      testthat::expect_equal(rr1, rr2)
    }
  })


testthat::test_that(
  "residuals_phregr: counting process data with covariates", {
    pbcseq <- pbcseq %>%
      group_by(id) %>%
      arrange(id, day) %>%
      mutate(tstart = day,
             tstop = ifelse(row_number() != n(), lead(day), futime),
             event = ifelse(row_number() != n(), 0, 1*(status == 2)))

    for (type in c("martingale", "deviance", "score", "schoenfeld",
                   "dfbeta", "dfbetas", "scaledsch")) {
      fit1 <- phregr(pbcseq, time="tstart", time2="tstop", event="event",
                     covariates=c("age", "edema", "log(bili)",
                                  "log(protime)", "log(albumin)"),
                     id="id", robust=TRUE)

      rr1 <- residuals_phregr(fit1, type=type)

      fit2 <- coxph(Surv(tstart, tstop, event) ~ age + edema + log(bili) +
                      log(protime) + log(albumin), data=pbcseq,
                    cluster=id, robust=TRUE)

      rr2 <- resid(fit2, type=type)

      if (type=="martingale" || type=="deviance") {
        names(rr2) <- NULL
      } else {
        rownames(rr2) <- NULL
        if (type=="schoenfeld" || type=="scaledsch") {
          attr(rr1, "time") <- NULL
        } else if (type=="dfbeta" || type=="dfbetas") {
          colnames(rr1) <- NULL
        }
      }

      testthat::expect_equal(rr1, rr2)
    }
  })

