#' @title Parametric Regression Models for Failure Time Data
#' @description Obtains the parameter estimates from parametric
#' regression models with uncensored, right censored, left censored, or
#' interval censored data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for interval censored data.
#'
#'   * \code{time2}: The right end of each interval for interval
#'     censored data.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates.
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID to group the score residuals
#'     in computing the robust sandwich variance.
#'
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for interval censored data in the input data.
#' @param time2 The name of the right end of each interval for
#'   interval censored data in the input data.
#' @param event The name of the event variable in the input data
#'   for right censored data.
#' @param covariates The vector of names of baseline covariates
#'   in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param dist The assumed distribution for time to event. Options include
#'   "exponential", "weibull", "lognormal", and "loglogistic" to be
#'   modeled on the log-scale, and "normal" and "logistic" to be modeled
#'   on the original scale.
#' @param init A vector of initial values for the model parameters,
#'   including regression coefficients and the log scale parameter.
#'   By default, initial values are derived from an intercept-only model.
#'   If this approach fails, ordinary least squares (OLS) estimates,
#'   ignoring censoring, are used instead.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence.
#'
#' @details There are two ways to specify the model, one for right censored
#' data through the time and event variables, and the other for interval
#' censored data through the time (lower) and time2 (upper) variables.
#' For the second form, we follow the convention used in SAS PROC LIFEREG:
#'
#' * If lower is not missing, upper is not missing, and lower is equal
#'   to upper, then there is no censoring and the event occurred at
#'   time lower.
#'
#' * If lower is not missing, upper is not missing, and lower < upper,
#'   then the event time is censored within the interval (lower, upper).
#'
#' * If lower is missing, but upper is not missing, then upper will be
#'   used as the left censoring value.
#'
#' * If lower is not missing, but upper is missing, then lower will be
#'   used as the right censoring value.
#'
#' * If lower is not missing, upper is not missing, but lower > upper,
#'   or if both lower and upper are missing, then the observation will
#'   not be used.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of observations.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum log-likelihood.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{dist}: The assumed distribution.
#'
#'     - \code{p}: The number of parameters, including the intercept,
#'       regression coefficients associated with the covariates, and
#'       the log scale parameters for the strata.
#'
#'     - \code{nvar}: The number of regression coefficients associated
#'       with the covariates (excluding the intercept).
#'
#'     - \code{robust}: Whether the robust sandwich variance estimate
#'       is requested.
#'
#'     - \code{fail}: Whether the model fails to converge.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The parameter estimate.
#'
#'     - \code{sebeta}: The standard error of parameter estimate.
#'
#'     - \code{z}: The Wald test statistic for the parameter.
#'
#'     - \code{expbeta}: The exponentiated parameter estimate.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of parameter estimate
#'       if robust variance is requested.
#'
#' * \code{linear_predictors}: The vector of linear predictors.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{nvar}: The number of columns of the design matrix excluding
#'   the intercept.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' John D. Kalbfleisch and Ross L. Prentice.
#' The Statistical Analysis of Failure Time Data.
#' Wiley: New York, 1980.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # right censored data
#' (fit1 <- liferegr(
#'   data = rawdata %>% filter(iterationNumber == 1) %>%
#'          mutate(treat = (treatmentGroup == 1)),
#'   stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", dist = "weibull"))
#'
#' # tobit regression for left censored data
#' (fit2 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal"))
#'
#' @export
liferegr <- function(data, stratum = "", time = "time", time2 = "",
                     event = "event", covariates = "", weight = "",
                     offset = "", id = "", dist = "weibull",
                     init = NA_real_, robust = FALSE, plci = FALSE,
                     alpha = 0.05, maxiter = 50, eps = 1.0e-9) {

  # validate input
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data frame");
  }

  if (inherits(data, "data.table") || inherits(data, "tbl") ||
      inherits(data, "tbl_df")) {
    df <- as.data.frame(data)
  } else {
    df <- data
  }

  for (nm in c(time, time2, event, weight, offset, id)) {
    if (!is.character(nm) || length(nm) != 1) {
      stop(paste(nm, "must be a single character string."));
    }
  }

  # select complete cases for the relevant variables
  elements <- unique(c(stratum, covariates, weight, offset, id))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)

  # check if the input data contains the required columns
  missing_cols <- setdiff(var_all, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("The following required columns are missing in the input data: ",
                paste(missing_cols, collapse = ", ")))
  }

  # use complete.cases on the subset of columns we care about
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  # Determine if covariates were provided (empty string or NULL means no covariates)
  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  # build design matrix and extract variable names
  if (misscovariates) {
    t1 <- terms(formula("~1"))
    param <- "(Intercept)"
    varnames <- ""
    xlevels <- NULL
  } else {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      param <- c("(Intercept)", covariates)
      varnames <- covariates
      t1 <- terms(fml_cov)
      xlevels <- NULL
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      param <- colnames(mm)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      t1 <- terms(fml_cov)
      xlevels <- mf$xlev
      # copy model-matrix columns into df only if they are missing
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  }

  # call the core fitting function
  fit <- liferegRcpp(
    data = df,
    stratum = stratum,
    time = time,
    time2 = time2,
    event = event,
    covariates = varnames,
    weight = weight,
    offset = offset,
    id = id,
    dist = dist,
    init = init,
    robust = robust,
    plci = plci,
    alpha = alpha,
    maxiter = maxiter,
    eps = eps)

  # post-process the output
  fit$p <- fit$sumstat$p[1]
  fit$nvar <- fit$sumstat$nvar[1]

  if (fit$p > 0) {
    par <- fit$parest$param[1:fit$p]
    if (length(par) > length(param)) {
      fit$param <- c(param, par[(1+length(param)):length(par)])
    } else {
      fit$param <- param
    }

    fit$beta <- fit$parest$beta
    names(fit$beta) <- fit$param

    dimnames(fit$vbeta) <- list(fit$param, fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) <- list(fit$param, fit$param)
    }
  }

  fit$terms <- t1
  if (fit$p > 0) fit$xlevels <- xlevels

  fit$settings <- list(
    data = data,
    stratum = stratum,
    time = time,
    time2 = time2,
    event = event,
    covariates = covariates,
    weight = weight,
    offset = offset,
    id = id,
    dist = dist,
    init = init,
    robust = robust,
    plci = plci,
    alpha = alpha,
    maxiter = maxiter,
    eps = eps
  )

  class(fit) <- "liferegr"
  fit
}



#' @title Proportional Hazards Regression Models
#' @description Obtains the hazard ratio estimates from the proportional
#' hazards regression model with right censored or counting process data.
#'
#' @param data The input data frame that contains the following variables:
#'
#'   * \code{stratum}: The stratum.
#'
#'   * \code{time}: The follow-up time for right censored data, or
#'     the left end of each interval for counting process data.
#'
#'   * \code{time2}: The right end of each interval for counting process
#'     data. Intervals are assumed to be open on the left
#'     and closed on the right, and event indicates whether an event
#'     occurred at the right end of each interval.
#'
#'   * \code{event}: The event indicator, 1=event, 0=no event.
#'
#'   * \code{covariates}: The values of baseline covariates (and
#'     time-dependent covariates in each interval for counting
#'     process data).
#'
#'   * \code{weight}: The weight for each observation.
#'
#'   * \code{offset}: The offset for each observation.
#'
#'   * \code{id}: The optional subject ID for counting process data
#'     with time-dependent covariates.
#'
#' @param stratum The name(s) of the stratum variable(s) in the input data.
#' @param time The name of the time variable or the left end of each
#'   interval for counting process data in the input data.
#' @param time2 The name of the right end of each interval for counting
#'   process data in the input data.
#' @param event The name of the event variable in the input data.
#' @param covariates The vector of names of baseline and time-dependent
#'   covariates in the input data.
#' @param weight The name of the weight variable in the input data.
#' @param offset The name of the offset variable in the input data.
#' @param id The name of the id variable in the input data.
#' @param ties The method for handling ties, either "breslow" or
#'   "efron" (default).
#' @param init The vector of initial values. Defaults to zero for all
#'   variables.
#' @param robust Whether a robust sandwich variance estimate should be
#'   computed. In the presence of the id variable, the score residuals
#'   will be aggregated for each id when computing the robust sandwich
#'   variance estimate.
#' @param est_basehaz Whether to estimate the baseline hazards.
#'   Defaults to \code{TRUE}.
#' @param est_resid Whether to estimate the martingale residuals.
#'   Defaults to \code{TRUE}.
#' @param firth Whether to use Firth’s penalized likelihood method.
#'   Defaults to \code{FALSE}.
#' @param plci Whether to obtain profile likelihood confidence interval.
#' @param alpha The two-sided significance level.
#' @param maxiter The maximum number of iterations.
#' @param eps The tolerance to declare convergence.
#'
#' @return A list with the following components:
#'
#' * \code{sumstat}: The data frame of summary statistics of model fit
#'   with the following variables:
#'
#'     - \code{n}: The number of observations.
#'
#'     - \code{nevents}: The number of events.
#'
#'     - \code{loglik0}: The (penalized) log-likelihood under null.
#'
#'     - \code{loglik1}: The maximum (penalized) log-likelihood.
#'
#'     - \code{scoretest}: The score test statistic.
#'
#'     - \code{niter}: The number of Newton-Raphson iterations.
#'
#'     - \code{ties}: The method for handling ties, either "breslow" or
#'       "efron".
#'
#'     - \code{p}: The number of columns of the Cox model design matrix.
#'
#'     - \code{robust}: Whether to use the robust variance estimate.
#'
#'     - \code{firth}: Whether to use Firth's penalized likelihood method.
#'
#'     - \code{fail}: Whether the model fails to converge.
#'
#'     - \code{loglik0_unpenalized}: The unpenalized log-likelihood under null.
#'
#'     - \code{loglik1_unpenalized}: The maximum unpenalized log-likelihood.
#'
#' * \code{parest}: The data frame of parameter estimates with the
#'   following variables:
#'
#'     - \code{param}: The name of the covariate for the parameter estimate.
#'
#'     - \code{beta}: The log hazard ratio estimate.
#'
#'     - \code{sebeta}: The standard error of log hazard ratio estimate.
#'
#'     - \code{z}: The Wald test statistic for log hazard ratio.
#'
#'     - \code{expbeta}: The hazard ratio estimate.
#'
#'     - \code{lower}: The lower limit of confidence interval.
#'
#'     - \code{upper}: The upper limit of confidence interval.
#'
#'     - \code{p}: The p-value from the chi-square test.
#'
#'     - \code{method}: The method to compute the confidence interval and
#'       p-value.
#'
#'     - \code{sebeta_naive}: The naive standard error of log hazard ratio
#'       estimate if robust variance is requested.
#'
#' * \code{basehaz}: The data frame of baseline hazards with the following
#'   variables (if est_basehaz is TRUE):
#'
#'     - \code{time}: The observed event time.
#'
#'     - \code{nrisk}: The number of patients at risk at the time point.
#'
#'     - \code{nevent}: The number of events at the time point.
#'
#'     - \code{haz}: The baseline hazard at the time point.
#'
#'     - \code{varhaz}: The variance of the baseline hazard at the time point
#'       assuming the parameter beta is known.
#'
#'     - \code{gradhaz}: The gradient of the baseline hazard with respect to
#'       beta at the time point (in the presence of covariates).
#'
#'     - \code{stratum}: The stratum.
#'
#' * \code{residuals}: The martingale residuals.
#'
#' * \code{linear_predictors}: The vector of linear predictors.
#'
#' * \code{p}: The number of parameters.
#'
#' * \code{param}: The parameter names.
#'
#' * \code{beta}: The parameter estimate.
#'
#' * \code{vbeta}: The covariance matrix for parameter estimates.
#'
#' * \code{vbeta_naive}: The naive covariance matrix for parameter estimates.
#'
#' * \code{terms}: The terms object.
#'
#' * \code{xlevels}: A record of the levels of the factors used in fitting.
#'
#' * \code{settings}: A list containing the input parameter values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Per K. Anderson and Richard D. Gill.
#' Cox's regression model for counting processes, a large sample study.
#' Annals of Statistics 1982; 10:1100-1120.
#'
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' (fit1 <- phregr(
#'   data = rawdata %>% filter(iterationNumber == 1) %>%
#'     mutate(treat = 1*(treatmentGroup == 1)),
#'   stratum = "stratum",
#'   time = "timeUnderObservation", event = "event",
#'   covariates = "treat", est_basehaz = FALSE, est_resid = FALSE))
#'
#' # Example 2 with counting process data and robust variance estimate
#' (fit2 <- phregr(
#'   data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'   time = "start", time2 = "stop", event = "event",
#'   covariates = c("rx", "age"), id = "id",
#'   robust = TRUE, est_basehaz = TRUE, est_resid = TRUE))
#'
#' @export
phregr <- function(data, stratum = "", time = "time", time2 = "",
                   event = "event", covariates = "", weight = "",
                   offset = "", id = "", ties = "efron",
                   init = NA_real_,  robust = FALSE,
                   est_basehaz = TRUE, est_resid = TRUE,
                   firth = FALSE, plci = FALSE, alpha = 0.05,
                   maxiter = 50, eps = 1.0e-9) {

  # validate input
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data frame");
  }

  if (inherits(data, "data.table") || inherits(data, "tbl") ||
      inherits(data, "tbl_df")) {
    df <- as.data.frame(data)
  } else {
    df <- data
  }

  for (nm in c(time, time2, event, weight, offset, id)) {
    if (!is.character(nm) || length(nm) != 1) {
      stop(paste(nm, "must be a single character string."));
    }
  }

  # select complete cases for the relevant variables
  elements <- unique(c(stratum, time, time2, event, covariates, weight, offset, id))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)

  # check if the input data contains the required columns
  missing_cols <- setdiff(var_all, names(df))
  if (length(missing_cols) > 0) {
    stop(paste0("The following required columns are missing in the input data: ",
                paste(missing_cols, collapse = ", ")))
  }

  # use complete.cases on the subset of columns we care about
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  # Determine if covariates were provided (empty string or NULL means no covariates)
  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  # build design matrix and extract variable names
  if (misscovariates) {
    t1 <- terms(formula("~1"))
    param <- "(Intercept)"
    varnames <- ""
    xlevels <- NULL
  } else {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      param <- c("(Intercept)", covariates)
      varnames <- covariates
      t1 <- terms(fml_cov)
      xlevels <- NULL
    } else {
      # Use model.matrix to handle factors and interactions
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      param <- colnames(mm)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      t1 <- terms(fml_cov)
      xlevels <- mf$xlev
      # Add derived columns to df if not already present
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  }

  # call the core fitting function
  fit <- phregRcpp(
    data = df,
    stratum = stratum,
    time = time,
    time2 = time2,
    event = event,
    covariates = varnames,
    weight = weight,
    offset = offset,
    id = id,
    ties = ties,
    init = init,
    robust = robust,
    est_basehaz = est_basehaz,
    est_resid = est_resid,
    firth = firth,
    plci = plci,
    alpha = alpha,
    maxiter = maxiter,
    eps = eps)

  # post-process the output
  fit$p <- fit$sumstat$p[1]

  if (fit$p > 0) {
    fit$param <- param[-1]
    fit$beta <- fit$parest$beta
    names(fit$beta) <- fit$param

    dimnames(fit$vbeta) <- list(fit$param, fit$param)
    if (robust) {
      dimnames(fit$vbeta_naive) <- list(fit$param, fit$param)
    }
  }

  fit$terms <- t1
  if (fit$p > 0) fit$xlevels <- xlevels

  fit$settings <- list(
    data = data,
    stratum = stratum,
    time = time,
    time2 = time2,
    event = event,
    covariates = covariates,
    weight = weight,
    offset = offset,
    id = id,
    ties = ties,
    iniy = init,
    robust = robust,
    est_basehaz = est_basehaz,
    est_resid = est_resid,
    firth = firth,
    plci = plci,
    alpha = alpha,
    maxiter = maxiter,
    eps = eps
  )

  class(fit) <- "phregr"
  fit
}



#' @title Survival Curve for Proportional Hazards Regression Models
#' @description Obtains the predicted survivor function for a proportional
#' hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
#' @param newdata A data frame with the same variable names as those that
#'   appear in the \code{phregr} call. For right-censored data, one curve is
#'   produced per row to represent a cohort whose covariates correspond to
#'   the values in \code{newdata}. For counting-process data, one curve is
#'   produced per \code{id} in \code{newdata} to present the survival curve
#'   along the path of time-dependent covariates at the observed event
#'   times in the data used to fit \code{phregr}.
#' @param sefit Whether to compute the standard error of the survival
#'   estimates.
#' @param conftype The type of the confidence interval. One of \code{"none"},
#'   \code{"plain"}, \code{"log"}, \code{"log-log"} (the default), or
#'   \code{"arcsin"}. The \code{arcsin} option bases the intervals on
#'   \code{asin(sqrt(surv))}.
#' @param conflev The level of the two-sided confidence interval for
#'   the survival probabilities. Defaults to 0.95.
#'
#' @details
#' If \code{newdata} is not provided and there is no covariate, survival
#' curves based on the \code{basehaz} data frame will be produced.
#'
#' @return A data frame with the following variables:
#'
#' * \code{id}: The id of the subject for counting-process data with
#'   time-dependent covariates.
#'
#' * \code{time}: The observed times in the data used to fit
#'   \code{phregr}.
#'
#' * \code{nrisk}: The number of patients at risk at the time point in the
#'   data used to fit \code{phregr}.
#'
#' * \code{nevent}: The number of patients having event at the time point
#'   in the data used to fit \code{phregr}.
#'
#' * \code{cumhaz}: The cumulative hazard at the time point.
#'
#' * \code{surv}: The estimated survival probability at the time point.
#'
#' * \code{sesurv}: The standard error of the estimated survival probability.
#'
#' * \code{lower}: The lower confidence limit for survival probability.
#'
#' * \code{upper}: The upper confidence limit for survival probability.
#'
#' * \code{conflev}: The level of the two-sided confidence interval.
#'
#' * \code{conftype}: The type of the confidence interval.
#'
#' * \code{covariates}: The values of covariates based on \code{newdata}.
#'
#' * \code{stratum}: The stratum of the subject.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau and Patricia M. Grambsch.
#' Modeling Survival Data: Extending the Cox Model.
#' Springer-Verlag, 2000.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
#'                  mutate(treat = 1*(treatmentGroup == 1)),
#'                stratum = "stratum",
#'                time = "timeUnderObservation", event = "event",
#'                covariates = "treat")
#'
#' surv1 <- survfit_phregr(fit1,
#'                         newdata = data.frame(
#'                           stratum = as.integer(c(1,1,2,2)),
#'                           treat = c(1,0,1,0)))
#' head(surv1)
#'
#' # Example 2 with counting process data and robust variance estimate
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' surv2 <- survfit_phregr(fit2,
#'                         newdata = data.frame(
#'                           id = c(4,4,11,11),
#'                           age = c(-7.737,-7.737,-0.019,-0.019),
#'                           start = c(0,36,0,26),
#'                           stop = c(36,39,26,153),
#'                           rx = c(0,1,0,1)))
#' head(surv2)
#'
#' @export
survfit_phregr <- function(object, newdata, sefit = TRUE,
                           conftype = "log-log", conflev = 0.95) {

  if (!inherits(object, "phregr")) stop("object must be of class 'phregr'");

  p <- object$p
  if (p == 0) {
    beta <- 0.0
    vbeta <- 0.0
  } else {
    beta <- object$beta
    vbeta <- object$vbeta
  }

  basehaz <- object$basehaz

  covariates <- object$settings$covariates
  stratum <- object$settings$stratum
  offset <- object$settings$offset
  id <- object$settings$id

  if (id != "") {
    tstart <- object$settings$time
    tstop <- object$settings$time2
  } else {
    tstart <- ""
    tstop <- ""
  }

  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  if (misscovariates && !(missing(newdata) || is.null(newdata))) {
    stop("covariates must be specified when newdata is available")
  }

  if (!(missing(newdata) || is.null(newdata))) {
    df <- newdata
    fml_cov <- formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- covariates
    } else {
      # Use model.matrix to handle factors and interactions
      mf <- model.frame(fml_cov, data = df, na.action = na.pass, xlev = object$xlevels)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    if (p > 0) {
      stop("newdata must be provided for Cox models with covariates")
    } else {
      p_stratum <- length(stratum);
      if (p_stratum == 0 || (p_stratum == 1 && stratum[1] == "")) {
        df <- data.frame(dummy_x_ = 0)
      } else {
        df <- unique(basehaz[, stratum, drop = FALSE])
      }
    }

    beta <- NA
    vbeta <- NA
    varnames <- ""
  }

  if (!is.matrix(vbeta)) vbeta <- as.matrix(vbeta)

  if (missing(basehaz) || is.null(basehaz)) {
    stop("basehaz must be provided")
  }

  survfit_phregRcpp(p = p, beta = beta, vbeta = vbeta, basehaz = basehaz,
                    newdata = df, covariates = varnames,
                    stratum = stratum, offset = offset, id = id,
                    tstart = tstart, tstop = tstop, sefit = sefit,
                    conftype = conftype, conflev = conflev)
}



#' @title Residuals for Parametric Regression Models for Failure Time Data
#' @description Obtains the response, martingale, deviance, dfbeta, and
#' likelihood displacement residuals for a parametric regression model
#' for failure time data.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"response"}, \code{"martingale"}, \code{"deviance"},
#'   \code{"dfbeta"}, \code{"dfbetas"}, \code{"working"}, \code{"ldcase"},
#'   \code{"ldresp"}, \code{"ldshape"}, and \code{"matrix"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' The algorithms follow the \code{residuals.survreg} function in the
#' \code{survival} package, except for martingale residuals, which
#' are defined only for event or right-censored data for exponential,
#' weibull, lognormal, and loglogistic distributions.
#'
#' @return
#' Either a vector or a matrix of residuals, depending on the specified type:
#'
#' * \code{response} residuals are on the scale of the original data.
#'
#' * \code{martingale} residuals are event indicators minus the cumulative
#'   hazards for event or right-censored data.
#'
#' * \code{working} residuals are on the scale of the linear predictor.
#'
#' * \code{deviance} residuals are on the log-likelihood scale.
#'
#' * \code{dfbeta} residuals are returned as a matrix, where the
#'   \eqn{i}-th row represents the approximate change in the model
#'   coefficients resulting from the inclusion of subject \eqn{i}.
#'
#' * \code{dfbetas} residuals are similar to \code{dfbeta} residuals, but
#'   each column is scaled by the standard deviation of the
#'   corresponding coefficient.
#'
#' * \code{matrix} residuals are a matrix of derivatives of the
#'   log-likelihood function. Let \eqn{L} be the log-likelihood, \eqn{p} be
#'   the linear predictor (\eqn{X\beta}), and \eqn{s} be \eqn{log(\sigma)}.
#'   Then the resulting matrix contains six columns: \eqn{L},
#'   \eqn{\partial L/\partial p}, \eqn{\partial^2 L/\partial p^2},
#'   \eqn{\partial L/\partial s}, \eqn{\partial^2 L/\partial s^2}, and
#'   \eqn{\partial L^2/\partial p\partial s}.
#'
#' * \code{ldcase} residulas are likelihood displacement for case weight
#'   perturbation.
#'
#' * \code{ldresp} residuals are likelihood displacement for response value
#'   perturbation.
#'
#' * \code{ldshape} residuals are likelihood displacement related to the
#'   shape parameter.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Escobar, L. A. and Meeker, W. Q.
#' Assessing influence in regression analysis with censored data.
#' Biometrics 1992; 48:507-528.
#'
#' @examples
#'
#' library(dplyr)
#'
#' fit1 <- liferegr(
#'   data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
#'   time = "time", time2 = "durable",
#'   covariates = c("age", "quant"), dist = "normal")
#'
#' resid <- residuals_liferegr(fit1, type = "response")
#' head(resid)
#'
#' @export
residuals_liferegr <- function(
    object, type = c("response", "martingale", "deviance", "dfbeta", "dfbetas",
                     "working", "ldcase", "ldresp", "ldshape", "matrix"),
    collapse = FALSE, weighted = (type %in% c("dfbeta", "dfbetas"))) {

  if (!inherits(object, "liferegr")) stop("object must be of class 'liferegr'");

  p <- object$p
  df <- object$settings$data
  stratum <- object$settings$stratum
  covariates <- object$settings$covariates
  weight <- object$settings$weight
  offset <- object$settings$offset
  id <- object$settings$id

  elements <- unique(c(stratum, covariates, weight, offset, id))
  elements <- elements[elements != ""]
  if (!(length(elements) == 0)) {
    fml_all <- formula(paste("~", paste(elements, collapse = "+")))
    var_all <- all.vars(fml_all)
    rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
    if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
    df <- df[rows_ok, , drop = FALSE]
  }

  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  if (!misscovariates) {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- covariates
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }

  type <- match.arg(type)

  if (type == "dfbeta" || type == "dfbetas") {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }

  vv <- drop(object$vbeta_naive)
  if (is.null(vv)) vv <- drop(object$vbeta)

  rr <- residuals_liferegRcpp(beta = object$beta,
                              vbeta = vv,
                              data = df,
                              stratum = stratum,
                              time = object$settings$time,
                              time2 = object$settings$time2,
                              event = object$settings$event,
                              covariates = varnames,
                              weight = weight,
                              offset = offset,
                              id = id,
                              dist = object$settings$dist,
                              type = type,
                              collapse = collapse,
                              weighted = weighted)

  if (type == "response" || type == "martingale" || type == "deviance" ||
      type == "working" || type == "ldcase" || type == "ldresp" ||
      type == "ldshape") {
    rr <- as.numeric(rr)
  } else if (type == "dfbeta" || type == "dfbetas") {
    colnames(rr) <- object$param
  } else {
    colnames(rr) <- c("g", "dg", "ddg", "ds", "dds", "dsg");
  }

  rr
}


#' @title Residuals for Proportional Hazards Regression Models
#' @description Obtains the martingale, deviance, score, or Schoenfeld
#' residuals for a proportional hazards regression model.
#'
#' @param object The output from the \code{phregr} call.
#' @param type The type of residuals desired, with options including
#'   \code{"martingale"}, \code{"deviance"}, \code{"score"},
#'   \code{"schoenfeld"}, \code{"dfbeta"}, \code{"dfbetas"}, and
#'   \code{"scaledsch"}.
#' @param collapse Whether to collapse the residuals by \code{id}.
#'   This is not applicable for Schoenfeld type residuals.
#' @param weighted Whether to compute weighted residuals.
#'
#' @details
#' For score and Schoenfeld type residuals, the proportional hazards model
#' must include at least one covariate. The algorithms for \code{deviance},
#' \code{dfbeta}, \code{dfbetas}, and \code{scaledsch} residuals follow
#' the \code{residuals.coxph} function in the \code{survival} package.
#'
#' @return For martingale and deviance residuals, the result is a vector
#' with one element corresponding to each subject (without \code{collapse}).
#' For score residuals, the result is a matrix where each row represents
#' a subject and each column corresponds to a variable. The row order
#' aligns with the input data used in the original fit. For Schoenfeld
#' residuals, the result is a matrix with one row for each event and
#' one column per variable. These rows are sorted by time within strata,
#' with the attributes \code{stratum} and \code{time} included.
#'
#' Score residuals represent each individual's contribution to the score
#' vector. Two commonly used transformations of this are \code{dfbeta},
#' which represents the approximate change in the coefficient vector
#' if the observation is excluded, and \code{dfbetas}, which gives the
#' approximate change in the coefficients scaled by the standard error
#' of the coefficients.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' Terry M. Therneau, Patricia M. Grambsch, and Thomas M. Fleming.
#' Martingale based residuals for survival models.
#' Biometrika 1990; 77:147-160.
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
#'
#' @examples
#'
#' library(dplyr)
#'
#' # Example 1 with right-censored data
#' fit1 <- phregr(data = rawdata %>% filter(iterationNumber == 1) %>%
#'                  mutate(treat = 1*(treatmentGroup == 1)),
#'                stratum = "stratum",
#'                time = "timeUnderObservation", event = "event",
#'                covariates = "treat")
#'
#' ressco <- residuals_phregr(fit1, type = "score")
#' head(ressco)
#'
#' # Example 2 with counting process data
#' fit2 <- phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
#'                time = "start", time2 = "stop", event = "event",
#'                covariates = c("rx", "age"), id = "id", robust = TRUE)
#'
#' resssch <- residuals_phregr(fit2, type = "scaledsch")
#' head(resssch)
#'
#' @export
residuals_phregr <- function(
    object, type = c("martingale", "deviance", "score", "schoenfeld",
                     "dfbeta", "dfbetas", "scaledsch"),
    collapse = FALSE, weighted = (type %in% c("dfbeta", "dfbetas"))) {

  if (!inherits(object, "phregr")) stop("object must be of class 'phregr'");

  p <- object$p
  beta <- object$beta
  df <- object$settings$data
  stratum <- object$settings$stratum
  time <- object$settings$time
  event <- object$settings$event
  covariates <- object$settings$covariates
  weight <- object$settings$weight
  offset <- object$settings$offset
  id <- object$settings$id

  elements <- unique(c(stratum, time, event, covariates, weight, offset, id))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  if (!misscovariates) {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- covariates
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }

  type <- match.arg(type)

  if (type == "dfbeta" || type == "dfbetas") {
    if (missing(weighted))
      weighted <- TRUE  # different default for this case
  }

  vv <- object$vbeta_naive
  if (is.null(vv)) vv <- object$vbeta

  if (p == 0) { # null Cox model case
    beta <- 0
    vv <- matrix(0,1,1)
  }

  temp <- residuals_phregRcpp(p = p,
                              beta = beta,
                              vbeta = vv,
                              resmart = object$residuals,
                              data = df,
                              stratum = stratum,
                              time = time,
                              time2 = object$settings$time2,
                              event = event,
                              covariates = varnames,
                              weight = weight,
                              offset = offset,
                              id = id,
                              ties = object$settings$ties,
                              type = type,
                              collapse = collapse,
                              weighted = weighted)

  if (type == "martingale" || type == "deviance") {
    rr <- temp$resid
  } else {
    if (p == 1) {
      rr <- c(temp$resid)
    } else {
      rr <- temp$resid
    }

    if (type == "schoenfeld" || type == "scaledsch") {
      attr(rr, "time") <- temp$time
      if (length(temp) == 3) {
        attr(rr, "strata") <- temp$strata
      }
    }
  }

  if (is.matrix(rr)) colnames(rr) <- object$param

  rr
}


#' @title Assess Proportional Hazards Assumption Based on Supremum Test
#' @description Obtains the standardized score processes and the simulated
#' distribution under the null hypothesis as well as the p-values for
#' the supremum tests.
#'
#' @param object The output from the \code{phregr} call.
#' @param resample The number of simulation samples for the supremem test.
#' @param seed The random seed for the simulations.
#'
#' @details
#' The supremum test corresponds to the ASSESS statement with \code{ph}
#' option of SAS PROC PHREG.
#'
#' @return A list with the following components:
#'
#' * \code{time} the unique event times.
#'
#' * \code{score_t} the observed standardized score process.
#'
#' * \code{score_t_list} a list of simulated standardized score processes
#'   under the null hypothesis.
#'
#' * \code{max_abs_value} the supremum of the absolute value of the observed
#'   standardized score process for each covariate and the supremum of
#'   the sum of absolute values of the observed standardized score processes
#'   across all covariates.
#'
#' * \code{p_value} the p-values for the supremum tests for each covariate
#'   and the global test.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#' D. Y. Lin, L. J. Wei, and Z. Ying.
#' Checking the Cox model with cumulative sums of martingale-based
#' residuals.
#' Biometrika 1993; 80:557-572.
#'
#' @examples
#'
#' fit <- phregr(data = liver, time = "Time", event = "Status",
#'               covariates = c("log(Bilirubin)", "log(Protime)",
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'
#' aph <- assess_phregr(fit, resample = 1000, seed = 314159)
#'
#' aph
#'
#' plot(aph, nsim = 20)
#'
#' @export
assess_phregr <- function(object, resample = 1000, seed = 12345) {

  if (!inherits(object, "phregr")) stop("object must be of class 'phregr'");

  p <- object$p
  df <- object$settings$data
  stratum <- object$settings$stratum
  time <- object$settings$time
  event <- object$settings$event
  covariates <- object$settings$covariates
  weight <- object$settings$weight
  offset <- object$settings$offset

  elements <- unique(c(stratum, time, event, covariates, weight, offset))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  if (!misscovariates) {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- covariates
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }

  aph <- assess_phregRcpp(p = p,
                          beta = object$beta,
                          vbeta = object$vbeta,
                          data = df,
                          stratum = stratum,
                          time = time,
                          time2 = object$settings$time2,
                          event = event,
                          covariates = varnames,
                          weight = weight,
                          offset = offset,
                          ties = object$settings$ties,
                          resample = resample,
                          seed = seed)

  aph$covariates <- object$param

  class(aph) <- "assess_phregr"
  aph
}


#' @title Assess Proportional Hazards Assumption Based on Scaled
#' Schoenfeld Residuals
#' @description Obtains the scaled Schoenfeld residuals and tests the
#' proportional hazards assumption using a score test for the interaction
#' between each covariate and a transformed time variable.
#'
#' @param object The output from the \code{phregr} call.
#' @param transform A character string indicating how survival times
#'   should be transformed before the test is performed. Supported values
#'   include "identity", "log", "rank", and "km" (default).
#'
#' @details
#' This corresponds to the \code{cox.zph} function from the \code{survival}
#' package with \code{terms = FALSE} and \code{global = TRUE}.
#'
#' @return A list with the following components:
#'
#' * \code{table} A matrix with one row for each parameter and a final
#'   row for the global test. The columns contain the score test
#'   for adding the time-dependent term, the degrees of freedom,
#'   and the two-sided p-value.
#'
#' * \code{x} The transformed time values.
#'
#' * \code{time} The original (untransformed) event times, with tied event
#'   times repeated.
#'
#' * \code{strata} The stratum index for each event.
#'
#' * \code{y} The matrix of scaled Schoenfeld residuals, with one column
#'   for each parameter and one row for each event. Column names correspond
#'   to the parameter names.
#'
#' * \code{var} An approximate covariance matrix of the scaled Schoenfeld
#'   residuals, used to construct an approximate standard error band for
#'   plots.
#'
#' * \code{transform} the transformation applied to the time values.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @references
#'
#' Patricia M. Grambsch and Terry M. Therneau.
#' Proportional hazards tests and diagnostics based on weighted residuals.
#' Biometrika 1994; 81:515-26.
#'
#' @examples
#'
#' fit <- phregr(data = liver, time = "Time", event = "Status",
#'               covariates = c("log(Bilirubin)", "log(Protime)",
#'                              "log(Albumin)", "Age", "Edema"),
#'               ties = "breslow")
#'
#' zph <- zph_phregr(fit, transform = "log")
#'
#' zph$table
#'
#' @export
zph_phregr <- function(object, transform = "km") {

  if (!inherits(object, "phregr")) stop("object must be of class 'phregr'");

  p <- object$p
  df <- object$settings$data
  stratum <- object$settings$stratum
  time <- object$settings$time
  event <- object$settings$event
  covariates <- object$settings$covariates
  weight <- object$settings$weight
  offset <- object$settings$offset

  elements <- unique(c(stratum, time, event, covariates, weight, offset))
  elements <- elements[elements != ""]
  fml_all <- formula(paste("~", paste(elements, collapse = "+")))
  var_all <- all.vars(fml_all)
  rows_ok <- which(complete.cases(df[, var_all, drop = FALSE]))
  if (length(rows_ok) == 0) stop("No complete cases found for the specified variables.")
  df <- df[rows_ok, , drop = FALSE]

  misscovariates <- length(covariates) == 0 ||
    (length(covariates) == 1 && (covariates[1] == ""))

  if (!misscovariates) {
    fml_cov <- as.formula(paste("~", paste(covariates, collapse = "+")))

    # QUICK PATH: if all covariates present in df and are numeric, avoid model.matrix
    cov_present <- covariates %in% names(df)
    all_numeric <- FALSE
    if (all(cov_present)) {
      all_numeric <- all(vapply(df[ covariates ], is.numeric, logical(1)))
    }

    if (all_numeric) {
      # Build design columns directly from numeric covariates (intercept + columns)
      # This avoids model.matrix and is valid when covariates are simple numeric columns.
      varnames <- covariates
    } else {
      # FALLBACK (existing robust behavior): use model.frame + model.matrix on df
      mf <- model.frame(fml_cov, data = df, na.action = na.pass)
      mm <- model.matrix(fml_cov, mf)
      colnames(mm) <- make.names(colnames(mm))
      varnames <- colnames(mm)[-1]
      missing_cols <- setdiff(varnames, names(df))
      if (length(missing_cols) > 0) {
        for (vn in missing_cols) df[[vn]] <- mm[, vn, drop = TRUE]
      }
    }
  } else {
    varnames <- ""
  }

  zph <- zph_phregRcpp(p = p,
                       beta = object$beta,
                       vbeta = object$vbeta,
                       resmart = object$residuals,
                       data = df,
                       stratum = stratum,
                       time = time,
                       time2 = object$settings$time2,
                       event = event,
                       covariates = varnames,
                       weight = weight,
                       offset = offset,
                       ties = object$settings$ties,
                       transform = transform)

  rownames(zph$table) <- c(object$param, "GLOBAL")
  colnames(zph$table) <- c("chisq", "df", "p")
  colnames(zph$y) <- object$param

  class(zph) <- "cox.zph"
  zph
}


