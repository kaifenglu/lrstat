# Log-Rank Test Simulation for Phase 2/3 Seamless Design

Simulate phase 2/3 seamless design using a weighted log-rank test.
Analyses can be triggered either by the cumulative number of events
(combined for an active arm and the common control) or by pre-specified
calendar times.

## Usage

``` r
lrsim_seamless(
  M = 2,
  K = 1,
  criticalValues = NA,
  futilityBounds = NULL,
  hazardRatioH0s = 1,
  allocations = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  lambdas = NULL,
  gammas = NULL,
  n = NA,
  followupTime = NA,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 0,
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

  Number of active treatment arms in Phase 2.

- K:

  Number of sequential looks in Phase 3.

- criticalValues:

  Numeric vector of length \\K + 1\\ giving the critical value for the
  Wald statistic at each look (Look 1 through Look \\K + 1\\). Decision
  rule:

  - At Look 1, compute the Wald statistic for each active arm versus the
    common control. If the maximum of these statistics exceeds the Look
    1 critical value, stop for efficacy.

  - If the Look 1 stopping rule is not met, select the active arm with
    the largest Wald statistic as the "best" arm and continue with that
    arm only versus control at subsequent looks.

  - For each look \\j = 2,\ldots,K+1\\, compare the selected arm to
    control; if its Wald statistic exceeds the Look \\j\\ critical
    value, stop for efficacy; otherwise continue.

  - If no critical value is exceeded by Look \\K + 1\\, the procedure
    ends without rejection.

- futilityBounds:

  Numeric vector of length \\K\\ giving the futility boundaries for
  Phase 2 and the first \\K-1\\ looks in Phase 3. The study stops for
  futility:

  - in Phase 2 if all active treatment arms cross the phase-2 futility
    boundary;

  - in Phase 3 if the selected arm crosses the futility boundary at an
    interim look; If omitted, no interim futility stopping is applied.

- hazardRatioH0s:

  Numeric vector of length \\M\\. Hazard ratios under \\H_0\\ for each
  active arm versus the common control. Defaults to 1 for superiority
  tests.

- allocations:

  Integer or integer vector of length \\M + 1\\. Number of subjects per
  arm within a randomization block. A single value implies equal
  allocation; defaults to 1.

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

  List of length \\M\\ (one element per arm). Each element is a scalar
  or a numeric vector of event hazard rates for the corresponding arm,
  given by analysis interval and stratum as required by the simulation.

- gammas:

  List of length \\M\\ (one element per arm). Each element is a scalar
  or a numeric vector of dropout hazard rates for the corresponding arm,
  by analysis interval and stratum.

- n:

  Planned total sample size across all active arms and control.

- followupTime:

  Follow-up time for the last enrolled subject.

- fixedFollowup:

  Whether a fixed follow-up design is used. Defaults to `FALSE` for
  variable follow-up.

- rho1:

  The first parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- rho2:

  The second parameter of the Fleming-Harrington family of weighted
  log-rank test. Defaults to 0 for conventional log-rank test.

- plannedEvents:

  Numeric vector of length \\K + 1\\ giving the planned cumulative
  number of events to trigger Look 1 through Look \\K + 1\\. Each entry
  refers to the combined events for the first active arm and the common
  control. Use `plannedEvents` to schedule event-driven looks.

- plannedTime:

  Numeric vector of calendar times for the analyses. If `plannedTime` is
  supplied, analyses are scheduled by calendar time and `plannedEvents`
  should be left missing.

- maxNumberOfIterations:

  Number of Monte Carlo replications. Defaults to 1000.

- maxNumberOfRawDatasetsPerStage:

  Number of subject-level raw datasets to retain per stage (for selected
  replications).

- seed:

  Random seed for reproducibility.

- nthreads:

  Number of threads for parallel simulation. Use 0 to accept the default
  RcppParallel behavior.

## Value

An S3 object of class `"lrsim_seamless"` with these components:

- `overview`: A list summarizing trial-level results and settings:

  - `selectAsBest`: Probability of selecting each active arm as the best
    arm at the end of phase 2.

  - `rejectPerStage`: Probability of rejecting the null for each active
    arm at each stage.

  - `futilityPerStage`: Probability of futility stopping for each active
    arm at each stage.

  - `cumulativeRejection`: Cumulative probability of rejection by stage.

  - `cumulativeFutility`: Cumulative futility stopping probabilities by
    stage.

  - `numberOfEvents`: Cumulative event counts by stage, including events
    from all arms in stage 1 and events from the selected arm and
    control in later stages.

  - `numberOfDropouts`: Cumulative dropouts by stage.

  - `numberOfSubjects`: Cumulative enrollments by stage.

  - `analysisTime`: Average calendar time for each stage among
    replications that reached that stage.

  - `overallReject`: Overall probability of rejecting the null by trial
    end.

  - `overallFutility`: Overall probability of stopping for futility by
    trial end.

  - `expectedNumberOfEvents`: Expected cumulative events at trial end.

  - `expectedNumberOfDropouts`: Expected cumulative dropouts at trial
    end.

  - `expectedNumberOfSubjects`: Expected cumulative enrollments at trial
    end.

  - `expectedStudyDuration`: Expected study duration.

  - `criticalValues`: The input critical values for each stage.

  - `futilityBounds`: The input futility boundaries for each stage.

  - `hazardRatioH0s`: The input hazard ratios under \\H_0\\.

  - `useEvents`: Logical indicating whether analyses were event-driven.

  - `numberOfIterations`: Number of simulation iterations performed.

  - `n`: Planned total sample size.

  - `fixedFollowup`: Logical indicating whether fixed follow-up was
    used.

  - `rho1`, `rho2`: Fleming–Harrington weighting parameters used.

  - `M`: Number of active arms in Phase 2.

  - `K`: Number of sequential looks in Phase 3.

- `sumdata1`: Data frame summarizing each iteration, stage, and
  treatment group:

  - `iterationNumber`, `eventsNotAchieved`, `stopStage`, `stageNumber`,
    `analysisTime`, `treatmentGroup`, `accruals`, `events`, `dropouts`.

  - For each stage the final row summarizes the overall study (all arms
    combined).

- `summdata2`: Data frame summarizing log-rank statistics by iteration,
  stage, and active arm:

  - `iterationNumber`, `bestArm`, `stopStage`, `stageNumber`,
    `analysisTime`, `activeArm`, `totalAccruals`, `totalEvents`,
    `totalDropouts`, `uscore`, `vscore`, `logRankStatistic`, `reject`,
    `futility`.

  - For each active arm, total accruals, events, and dropouts refer to
    the combined counts for that arm and the common control at that
    stage.

- `rawdata` (present when `maxNumberOfRawDatasetsPerStage` \> 0):
  Subject-level data for selected replications with variables:

  - `iterationNumber`, `stopStage`, `stageNumber`, `analysisTime`,
    `subjectId`, `arrivalTime`, `stratum`, `treatmentGroup`,
    `survivalTime`, `dropoutTime`, `timeUnderObservation`, `event`,
    `dropoutEvent`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(sim1 <- lrsim_seamless(
  M = 2,
  K = 2,
  criticalValues = c(3.852050, 2.723811, 2.223982),
  futilityBounds = c(0, 0.5),
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  lambdas = list(log(2)/12*0.5, log(2)/12*0.7, log(2)/12),
  n = 700,
  plannedEvents = c(40, 80, 120),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 0))
#> Error in print.lrsim_seamless(x): object 'M' not found
```
