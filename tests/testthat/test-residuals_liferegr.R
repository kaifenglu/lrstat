library(dplyr, warn.conflicts = FALSE)
library(survival)

testthat::test_that(
  "residuals_liferegr: right-censored data with covariates", {
    pbc <- pbc %>% mutate(event = 1*(status == 2))

    for (type in c("response", "deviance", "dfbeta", "dfbetas",
                   "working", "ldcase", "ldresp", "ldshape", "matrix")) {
      fit1 <- liferegr(pbc, time="time", event="event",
                       covariates=c("age", "edema", "log(bili)",
                                    "log(protime)", "log(albumin)"))

      rr1 <- residuals_liferegr(fit1, type=type)
      
      fit2 <- survreg(Surv(time, event) ~ age + edema + log(bili) +
                        log(protime) + log(albumin), data=pbc)

      rr2 <- resid(fit2, type=type)

      if (type=="response" || type=="deviance" || 
          type=="working" || type=="ldcase" ||
          type=="ldresp" || type=="ldshape") {
        names(rr2) <- NULL
      } else {
        rownames(rr2) <- NULL
        if (type=="dfbeta" || type=="dfbetas") {
          colnames(rr1) <- NULL
        }
      }

      testthat::expect_equal(rr1, rr2)
    }
  })
