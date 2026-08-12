# Log-Rank Test Simulation for MCPMod Design

Simulate an MCPMod design using a Cox model. Analyses can be triggered
either by the total number of events or by the pre-specified calendar
time.

## Usage

``` r
lrsim_mcpmod(
  M = 2,
  alpha = 0.05,
  hazardRatioH0s = 1,
  allocations = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  lambdas = NULL,
  candidateHazardRatios = NULL,
  gammas = NULL,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  plannedEvents = NA,
  plannedTime = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0
)
```

## Arguments

- M:

  Number of active treatment arms.

- alpha:

  Significance level for the max test (1-sided).

- hazardRatioH0s:

  Scalar or numeric vector of length \\M\\. Hazard ratios under \\H_0\\
  for each active arm versus the common control. Defaults to 1 for
  superiority tests.

- allocations:

  Integer or integer vector of length \\M + 1\\. Number of subjects per
  arm within a randomization block.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- piecewiseSurvivalTime:

  A vector that specifies the starting time of piecewise exponential
  survival time intervals. Must start with 0, e.g., `c(0, 6)` breaks the
  time axis into 2 event intervals: \\\[0, 6)\\ and \\\[6, \infty)\\.
  Defaults to 0 for exponential distribution.

- stratumFraction:

  A vector of stratum fractions that sum to 1. Defaults to 1 for no
  stratification.

- lambdas:

  List of length \\M + 1\\ (one element per arm). Each element is a
  scalar or a numeric vector of event hazard rates for the corresponding
  arm.

- candidateHazardRatios:

  Numeric matrix of dimension \\M \times T\\ containing the candidate
  hazard ratios for the \\T\\ candidate models.

- gammas:

  List of length \\M + 1\\ (one element per arm). Each element is a
  scalar or a numeric vector of dropout hazard rates for the
  corresponding arm.

- n:

  Planned total sample size across all active arms and control.

- followupTime:

  Follow-up time for the last enrolled subject.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- plannedEvents:

  Planned total number of events to trigger the analysis. Leave missing
  when using calendar-time planning.

- plannedTime:

  Planned calendar time for the analysis. Leave missing when using
  event-based planning.

- maxNumberOfIterations:

  Number of Monte Carlo replications.

- maxNumberOfRawDatasetsPerStage:

  Number of subject-level raw datasets to retain.

- seed:

  Random seed for reproducibility.

- nthreads:

  Number of threads for parallel simulation. Use 0 to accept the default
  RcppParallel behavior.

## Value

An S3 object of class `"lrsim_mcpmod"` with components `overview`,
`sumdata1`, `sumdata2`, and optionally `rawdata`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define candidate dose-response models
# (negative sign models decreasing dose-response for hazard)
f_emax <- function(d, ED50) { -d / (ED50 + d) }
f_exponential <- function(d, delta) { -exp(d / delta) + 1 }
f_linear <- function(d) { -d }
f_logistic <- function(d, ED50, delta) {
  -1 / (1 + exp((ED50 - d) / delta))
}
f_betamod <- function(d, delta1, delta2, D) {
  -dbeta(d / D, delta1 + 1, delta2 + 1)
}

# Log hazard parameters
b0 <- log(log(2) / 0.5)  # placebo: median survival 0.5 years
b1 <- log(0.6) + b0       # optimal dose: HR = 0.6 vs. placebo

# Rescale candidate models to match hazard rates
f1 <- function(dose) {
  a0 <- f_emax(0, 50); a1 <- f_emax(100, 50)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_emax(dose, 50)
}
f2 <- function(dose) {
  a0 <- f_exponential(0, 22.756); a1 <- f_exponential(100, 22.756)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_exponential(dose, 22.756)
}
f3 <- function(dose) {
  a0 <- f_emax(0, 6.25); a1 <- f_emax(100, 6.25)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_emax(dose, 6.25)
}
f4 <- function(dose) {
  a0 <- f_logistic(0, 40.3287, 6.9764)
  a1 <- f_logistic(100, 40.3287, 6.9764)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_logistic(dose, 40.3287, 6.9764)
}
f5 <- function(dose) {
  a0 <- f_linear(0); a1 <- f_linear(100)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_linear(dose)
}
f6 <- function(dose) {
  a0 <- f_betamod(0, 0.7489, 1.0485, 120)
  dstar <- 0.7489 * 120 / (0.7489 + 1.0485)
  a1 <- f_betamod(dstar, 0.7489, 1.0485, 120)
  slope <- (b1 - b0) / (a1 - a0)
  intercept <- b0 - slope * a0
  intercept + slope * f_betamod(dose, 0.7489, 1.0485, 120)
}

f <- list(f1, f2, f3, f4, f5, f6)
doselevels <- c(5, 25, 50, 100)

# Run the simulations (each candidate model is the true model in one scenario)
sims <- lapply(1:6, function(i) {
  lrsim_mcpmod(
    M = 4, alpha = 0.05, accrualIntensity = 360,
    lambdas = as.list(exp(f[[i]](c(doselevels, 0)))),
    candidateHazardRatios = exp(sapply(f, function(ff) ff(doselevels) - ff(0))),
    gammas = list(0, 0, 0, 0, 0),
    n = 300, plannedEvents = 242,
    maxNumberOfIterations = 1000,
    maxNumberOfRawDatasetsPerStage = 10,
    seed = 314159
  )
})

sapply(sims, function(s) s$overview$overallReject)
} # }
```
