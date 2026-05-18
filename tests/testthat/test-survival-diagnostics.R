testthat::test_that("diagnostics output object contracts and dimensions", {
  fit <- phregr(
    data = liver,
    time = "Time",
    event = "Status",
    covariates = c("log(Bilirubin)", "log(Protime)", "log(Albumin)", "Age", "Edema"),
    ties = "breslow"
  )

  zph <- zph_phregr(fit, transform = "log")
  testthat::expect_s3_class(zph, "cox.zph")
  testthat::expect_equal(colnames(zph$table), c("chisq", "df", "p"))
  testthat::expect_equal(nrow(zph$table), fit$p + 1L)
  testthat::expect_equal(ncol(zph$y), fit$p)
  testthat::expect_equal(dim(zph$var), c(fit$p, fit$p))
  testthat::expect_true(all(zph$table[, "p"] >= 0 & zph$table[, "p"] <= 1))

  aph <- assess_phregr(fit, resample = 30, seed = 123)
  testthat::expect_s3_class(aph, "assess_phregr")
  testthat::expect_named(
    aph,
    c("time", "score_t", "score_t_list", "max_abs_value", "p_value", "resample", "seed", "covariates")
  )
  testthat::expect_equal(ncol(aph$score_t), fit$p)
  testthat::expect_equal(length(aph$max_abs_value), fit$p + 1L)
  testthat::expect_equal(length(aph$p_value), fit$p + 1L)
  testthat::expect_true(all(aph$p_value >= 0 & aph$p_value <= 1))
  testthat::expect_equal(
    aph$max_abs_value,
    c(1.08798623329387, 1.72430587372516, 0.844321609048047,
      0.738742927033615, 1.43504821719989, 4.5107231623942),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    aph$p_value,
    c(0.0333333333333333, 0, 0.4, 0.566666666666667, 0.0333333333333333, 0),
    tolerance = 1e-12
  )

  sq <- survQuantile(
    time = c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6, 33.7, 17.1, 7.9, 10.5, 38),
    event = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1),
    probs = c(0.25, 0.5, 0.75)
  )
  testthat::expect_s3_class(sq, "data.frame")
  testthat::expect_named(sq, c("prob", "quantile", "lower", "upper", "cilevel", "transform"))
  testthat::expect_equal(nrow(sq), 3L)
  testthat::expect_true(all(sq$lower <= sq$quantile & sq$quantile <= sq$upper))
  testthat::expect_equal(sq$quantile, c(10.5, 38, 38), tolerance = 1e-12)
  testthat::expect_equal(sq$lower, c(3.9, 7.9, 19.5), tolerance = 1e-12)
  testthat::expect_equal(sq$upper, c(38, 38, 38), tolerance = 1e-12)

  kd <- kmdiff(
    data = rawdata[rawdata$iterationNumber == 1, ],
    stratum = "stratum",
    treat = "treatmentGroup",
    time = "timeUnderObservation",
    event = "event",
    milestone = 12
  )
  testthat::expect_s3_class(kd, "data.frame")
  testthat::expect_named(
    kd,
    c("milestone", "survDiffH0", "surv1", "surv2", "survDiff", "vsurv1", "vsurv2", "sesurvDiff",
      "survDiffZ", "survDiffPValue", "lower", "upper", "conflev")
  )
  testthat::expect_equal(nrow(kd), 1L)
  testthat::expect_equal(kd$survDiff, kd$surv1 - kd$surv2, tolerance = 1e-12)
  testthat::expect_equal(kd$surv1, 0.484842520079233, tolerance = 1e-12)
  testthat::expect_equal(kd$surv2, 0.385929654908856, tolerance = 1e-12)
  testthat::expect_equal(kd$survDiffZ, 2.20869280659367, tolerance = 1e-10)
  testthat::expect_equal(kd$survDiffPValue, 0.0271960151132257, tolerance = 1e-12)

  kmst <- kmstat(
    time = c(22, 40),
    milestone = 18,
    allocationRatioPlanned = 1,
    accrualTime = seq(0, 8),
    accrualIntensity = 26 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 6),
    stratumFraction = c(0.2, 0.8),
    lambda1 = c(0.0533, 0.0309, 1.5 * 0.0533, 1.5 * 0.0309),
    lambda2 = c(0.0533, 0.0533, 1.5 * 0.0533, 1.5 * 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 22,
    followupTime = 18,
    fixedFollowup = FALSE
  )
  testthat::expect_s3_class(kmst, "data.frame")
  testthat::expect_equal(nrow(kmst), 2L)
  testthat::expect_true(all(kmst$surv1 >= 0 & kmst$surv1 <= 1))
  testthat::expect_true(all(kmst$surv2 >= 0 & kmst$surv2 <= 1))
  testthat::expect_true(all(kmst$information > 0))
  testthat::expect_equal(kmst$subjects, c(468, 468), tolerance = 1e-12)
  testthat::expect_equal(kmst$surv1, c(0.384180499542573, 0.384180499542573), tolerance = 1e-12)
  testthat::expect_equal(kmst$surv2, c(0.266337409847818, 0.266337409847818), tolerance = 1e-12)
  testthat::expect_equal(kmst$survDiff, c(0.117843089694756, 0.117843089694756), tolerance = 1e-12)
  testthat::expect_equal(kmst$survDiffZ, c(1.33851229249307, 2.70570553501281), tolerance = 1e-10)

  ks <- kmsurv(
    time = c(2, 8),
    allocationRatioPlanned = 1,
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12
  )
  testthat::expect_equal(length(ks), 2L)
  testthat::expect_true(all(ks > 0 & ks < 1))
  testthat::expect_true(ks[1] > ks[2])
  testthat::expect_equal(ks, c(0.898885155151141, 0.667811626829863), tolerance = 1e-12)

  cr <- covrmst(
    t2 = 25,
    tau1 = 16,
    tau2 = 18,
    allocationRatioPlanned = 1,
    accrualTime = c(0, 3),
    accrualIntensity = c(10, 20),
    piecewiseSurvivalTime = c(0, 6),
    lambda1 = c(0.0533, 0.0309),
    lambda2 = c(0.0533, 0.0533),
    gamma1 = -log(1 - 0.05) / 12,
    gamma2 = -log(1 - 0.05) / 12,
    accrualDuration = 12,
    maxFollowupTime = 30
  )
  testthat::expect_equal(length(cr), 2L)
  testthat::expect_true(all(cr > 0))
  testthat::expect_equal(cr, c(0.373912675344442, 0.353265037559868), tolerance = 1e-12)

  u <- c(0, 1, 3, 4)
  lambda_pfs <- c(0.0151, 0.0403, 0.0501, 0.0558)
  lambda_os <- 0.0145
  hp <- hazard_pd(u, lambda_pfs, lambda_os, rho_pd_os = 0.5)
  testthat::expect_named(hp, c("piecewiseSurvivalTime", "hazard_pd", "hazard_os", "rho_pd_os"))
  testthat::expect_equal(length(hp$piecewiseSurvivalTime), length(hp$hazard_pd))
  testthat::expect_equal(length(hp$hazard_os), length(hp$hazard_pd))
  testthat::expect_true(all(hp$hazard_pd > 0))
  testthat::expect_equal(
    hp$piecewiseSurvivalTime,
    c(0, 1, 3, 3.19282466382887, 4, 5.38608514900017, 7.77912085911707,
      10.5416778452686, 13.8090892573467, 16.4219924831531, 17.8080776321533,
      22.9636703284218, 30.2300701153065, 42.6520625984596),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    hp$hazard_pd,
    c(0.00083567820519801, 0.0311408296162079, 0.0427651037242936,
      0.0430302690797641, 0.0498781128781187, 0.0505667622390376,
      0.0512705251354514, 0.051860161976481, 0.0522914905488956,
      0.0525291948431277, 0.0528407296754801, 0.0532889452117835,
      0.0537593531778292, 0.0541462469149736),
    tolerance = 1e-10
  )
  testthat::expect_equal(hp$hazard_os, rep(0.0145, 14), tolerance = 1e-12)
  testthat::expect_equal(hp$rho_pd_os, 0.5, tolerance = 1e-12)

  cp <- corr_pfs_os(u, lambda_pfs, lambda_os, rho_pd_os = 0.5)
  testthat::expect_equal(length(cp), 1L)
  testthat::expect_true(cp >= -1 && cp <= 1)
  testthat::expect_equal(cp, 0.548373565463488, tolerance = 1e-12)

  hs <- hazard_sub(
    u,
    hazard_itt = lambda_pfs,
    hazard_pos = c(0.0115, 0.0302, 0.0351, 0.0404),
    p_pos = 0.3
  )
  testthat::expect_named(hs, c("piecewiseSurvivalTime", "hazard_pos", "hazard_neg", "p_pos"))
  testthat::expect_equal(length(hs$piecewiseSurvivalTime), length(hs$hazard_neg))
  testthat::expect_equal(length(hs$hazard_pos), length(hs$hazard_neg))
  testthat::expect_true(all(hs$hazard_neg > 0))
  testthat::expect_equal(hs$piecewiseSurvivalTime, hp$piecewiseSurvivalTime, tolerance = 1e-10)
  testthat::expect_equal(
    hs$hazard_pos,
    c(0.0115, 0.0302, 0.0351, 0.0351, 0.0404, 0.0404, 0.0404,
      0.0404, 0.0404, 0.0404, 0.0404, 0.0404, 0.0404, 0.0404),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    hs$hazard_neg,
    c(0.0166468333509831, 0.0447145753828296, 0.0567660326689607,
      0.0568388172983284, 0.0628891862829182, 0.0631998134749536,
      0.0636527486677364, 0.0642314823290232, 0.0648500521955625,
      0.0653066964421317, 0.0661367543565381, 0.068011756180516,
      0.0722468751004948, 0.0816277208541779),
    tolerance = 1e-10
  )
  testthat::expect_equal(hs$p_pos, 0.3, tolerance = 1e-12)

  tcut <- c(0, 12, 36, 48)
  surv <- c(1, 0.95, 0.82, 0.74)
  lambda2 <- (log(surv[1:3]) - log(surv[2:4])) / (tcut[2:4] - tcut[1:3])

  bsim <- binary_tte_sim(
    kMax1 = 1,
    kMax2 = 2,
    accrualTime = seq(0, 8),
    accrualIntensity = 8 / 9 * seq(1, 9),
    piecewiseSurvivalTime = c(0, 12, 36),
    globalOddsRatio = 1,
    pi1 = 0.80,
    pi2 = 0.65,
    lambda1 = 0.65 * lambda2,
    lambda2 = lambda2,
    gamma1 = -log(1 - 0.04) / 12,
    gamma2 = -log(1 - 0.04) / 12,
    delta1 = -log(1 - 0.02) / 12,
    delta2 = -log(1 - 0.02) / 12,
    upper1 = 15 * 28 / 30.4,
    upper2 = 12 * 28 / 30.4,
    n = 80,
    plannedTime = 20 + 15 * 28 / 30.4,
    plannedEvents = c(20, 30),
    maxNumberOfIterations = 10,
    maxNumberOfRawDatasetsPerStage = 1,
    seed = 314159,
    nthreads = 1
  )
  testthat::expect_named(bsim, c("sumdataBIN", "sumdataTTE", "rawdataBIN", "rawdataTTE"))
  testthat::expect_s3_class(bsim$sumdataBIN, "data.frame")
  testthat::expect_s3_class(bsim$sumdataTTE, "data.frame")
  testthat::expect_true(nrow(bsim$sumdataBIN) > 0)
  testthat::expect_true(nrow(bsim$sumdataTTE) > 0)
  testthat::expect_true(all(c("riskDiff", "mhStatistic") %in% names(bsim$sumdataBIN)))
  testthat::expect_true(all(c("totalEvents", "logRankStatistic") %in% names(bsim$sumdataTTE)))
  testthat::expect_equal(dim(bsim$sumdataBIN), c(10L, 18L))
  testthat::expect_equal(dim(bsim$sumdataTTE), c(20L, 16L))
  testthat::expect_equal(dim(bsim$rawdataBIN), c(80L, 16L))
  testthat::expect_equal(dim(bsim$rawdataTTE), c(160L, 12L))
  testthat::expect_equal(bsim$sumdataBIN$analysisTime[1], 33.8157894736842, tolerance = 1e-10)
  testthat::expect_equal(bsim$sumdataBIN$riskDiff[1], 0.05, tolerance = 1e-12)
  testthat::expect_equal(bsim$sumdataBIN$seRiskDiff[1], 0.108108741552198, tolerance = 1e-12)
  testthat::expect_equal(bsim$sumdataBIN$mhStatistic[1], 0.462497290062881, tolerance = 1e-12)
  testthat::expect_equal(bsim$sumdataTTE$analysisTime[1], 68.3501167562583, tolerance = 1e-10)
  testthat::expect_equal(bsim$sumdataTTE$totalEvents[1], 20L)
  testthat::expect_equal(bsim$sumdataTTE$uscore[1], -2.22562303754909, tolerance = 1e-12)
  testthat::expect_equal(bsim$sumdataTTE$vscore[1], 4.98843571841856, tolerance = 1e-12)
  testthat::expect_equal(bsim$sumdataTTE$logRankStatistic[1], -0.996481907669256, tolerance = 1e-12)
})


testthat::test_that("selected outputs align with survival package references", {
  testthat::skip_if_not_installed("survival")

  fit <- phregr(
    data = liver,
    time = "Time",
    event = "Status",
    covariates = c("log(Bilirubin)", "log(Protime)", "log(Albumin)", "Age", "Edema"),
    ties = "breslow"
  )
  z_lr <- zph_phregr(fit, transform = "log")

  liver_ref <- transform(
    liver,
    logBilirubin = log(Bilirubin),
    logProtime = log(Protime),
    logAlbumin = log(Albumin)
  )
  fit_ref <- survival::coxph(
    survival::Surv(Time, Status) ~ logBilirubin + logProtime + logAlbumin + Age + Edema,
    data = liver_ref,
    ties = "breslow",
    x = TRUE,
    y = TRUE
  )
  z_ref <- survival::cox.zph(fit_ref, transform = "log", terms = FALSE, global = TRUE)

  testthat::expect_equal(unname(z_lr$table), unname(z_ref$table), tolerance = 1e-8)

  time_vec <- c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6, 33.7, 17.1, 7.9, 10.5, 38)
  event_vec <- c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1)
  probs <- c(0.25, 0.5, 0.75)

  sq <- survQuantile(time = time_vec, event = event_vec, probs = probs)

  sf <- survival::survfit(survival::Surv(time_vec, event_vec) ~ 1, conf.type = "log-log")
  q_ref <- as.numeric(quantile(sf, probs = probs)$quantile)
  testthat::expect_equal(sq$quantile, q_ref, tolerance = 1e-10)

  rd <- rawdata[rawdata$iterationNumber == 1, ]
  kd <- kmdiff(
    data = rd,
    stratum = "stratum",
    treat = "treatmentGroup",
    time = "timeUnderObservation",
    event = "event",
    milestone = 12
  )

  sf_group <- survival::survfit(
    survival::Surv(timeUnderObservation, event) ~ treatmentGroup,
    data = rd
  )
  s12 <- summary(sf_group, times = 12, extend = TRUE)$surv

  testthat::expect_equal(kd$surv1, s12[1], tolerance = 5e-4)
  testthat::expect_equal(kd$surv2, s12[2], tolerance = 5e-4)
  testthat::expect_equal(kd$survDiff, kd$surv1 - kd$surv2, tolerance = 1e-12)
})


testthat::test_that("non-PH synthetic data triggers expected diagnostic signals", {
  set.seed(20260518)

  n <- 600
  trt <- rbinom(n, 1, 0.5)
  tswitch <- 3

  # Crossing hazards: treatment has lower early hazard but higher late hazard.
  u <- runif(n)
  rate_early <- ifelse(trt == 1, 0.05, 0.18)
  rate_late <- ifelse(trt == 1, 0.25, 0.06)
  s_switch <- exp(-rate_early * tswitch)

  tevent <- ifelse(
    u > s_switch,
    -log(u) / rate_early,
    tswitch - log(u / s_switch) / rate_late
  )
  ctime <- rexp(n, rate = 0.01)
  time <- pmin(tevent, ctime)
  event <- as.integer(tevent <= ctime)

  d <- data.frame(time = time, event = event, trt = trt)

  fit_np <- phregr(
    data = d,
    time = "time",
    event = "event",
    covariates = "trt",
    ties = "breslow"
  )

  z_np <- zph_phregr(fit_np, transform = "log")
  testthat::expect_lt(z_np$table["trt", "p"], 1e-6)
  testthat::expect_lt(z_np$table["GLOBAL", "p"], 1e-6)

  aph_np <- assess_phregr(fit_np, resample = 120, seed = 99)
  testthat::expect_equal(length(aph_np$p_value), 2L)
  testthat::expect_lt(aph_np$p_value[1], 0.05)
  testthat::expect_lt(aph_np$p_value[2], 0.05)
})



testthat::test_that(
  "residuals_liferegr: right-censored data with covariates", {
    library(dplyr, warn.conflicts = FALSE)
    library(survival)
    pbc <- pbc %>% mutate(event = 1*(status == 2))

    fit1 <- liferegr(pbc, time="time", event="event",
                     covariates=c("age", "edema", "log(bili)",
                                  "log(protime)", "log(albumin)"))

    fit2 <- survreg(Surv(time, event) ~ age + edema + log(bili) +
                      log(protime) + log(albumin), data=pbc)

    for (type in c("response", "deviance", "dfbeta", "dfbetas",
                   "working", "ldcase", "ldresp", "ldshape", "matrix")) {

      rr1 <- residuals_liferegr(fit1, type=type)

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



testthat::test_that(
  "residuals_phregr: null model for right-censored data", {
    library(dplyr, warn.conflicts = FALSE)
    library(survival)
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
    library(dplyr, warn.conflicts = FALSE)
    library(survival)
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
    library(dplyr, warn.conflicts = FALSE)
    library(survival)
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
    library(dplyr, warn.conflicts = FALSE)
    library(survival)
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



testthat::test_that("survfit_phregr: right-censored data", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
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
           surv, sesurv, lower, upper)

  surv2b <- data.frame(time = surv2$time, nrisk = surv2$n.risk,
                       nevent = surv2$n.event, ncensor = surv2$n.censor,
                       cumhaz = surv2$cumhaz, surv = surv2$surv,
                       sesurv = surv2$surv*surv2$std.chaz,
                       lower = surv2$lower, upper = surv2$upper)

  testthat::expect_equal(surv1b, surv2b)
})


testthat::test_that("survfit_phregr: stratified analysis", {
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
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
  library(dplyr, warn.conflicts = FALSE)
  library(survival)
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
           surv, sesurv, lower, upper)

  surv2b <- data.frame(time = surv2$time, nrisk = surv2$n.risk,
                       nevent = surv2$n.event, ncensor = surv2$n.censor,
                       cumhaz = surv2$cumhaz, surv = surv2$surv,
                       sesurv = surv2$surv*surv2$std.chaz,
                       lower = surv2$lower, upper = surv2$upper)

  testthat::expect_equal(surv1b, surv2b)
})
