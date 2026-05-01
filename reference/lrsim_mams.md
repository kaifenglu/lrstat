# Log-Rank Test Simulation for Multi-Arm Multi-Stage Design

Simulate multi-arm multi-stage design using a weighted log-rank test.
Analyses can be triggered either by the cumulative number of events
(combined for an active arm and the common control) or by pre-specified
calendar times.

## Usage

``` r
lrsim_mams(
  M = 2,
  kMax = 1,
  criticalValues = NULL,
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

  Number of active treatment arms.

- kMax:

  Number of sequential looks.

- criticalValues:

  The matrix of by-level upper boundaries on the max z-test statistic
  scale for efficacy stopping. The first column is for level \\M\\, the
  second column is for level \\M - 1\\, and so on, with the last column
  for level 1. Decision rule:

  - At Look 1, compute the Wald statistic for each active arm versus the
    common control. If the maximum of these statistics exceeds the Look
    1 critical value, stop for efficacy and check whether there is any
    other active arm which can be rejected using a relaxed boundary
    under the closed testing principle.

  - If the Look 1 stopping rule is not met, continue to the next look;
    if the maximum of Wald statistics at this look exceeds the
    corresponding level \\M\\ critical value, stop for efficacy;
    otherwise continue.

  - If no critical value is exceeded by Look `kMax`, the procedure ends
    without rejection.

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

  - `rejectPerStage`: Probability of rejecting the null for each active
    arm at each stage.

  - `cumulativeRejection`: Cumulative probability of rejection for each
    active arm by stage.

  - `overallReject`: Overall probability of rejecting the null by trial
    end.

  - `numberOfEvents`: Cumulative event counts by stage and arm.

  - `numberOfDropouts`: Cumulative dropouts by stage and arm.

  - `numberOfSubjects`: Cumulative enrollments by stage and arm.

  - `analysisTime`: Average calendar time for each stage by arm among
    replications that reached that stage.

  - `expectedNumberOfEvents`: Expected cumulative events at trial end.

  - `expectedNumberOfDropouts`: Expected cumulative dropouts at trial
    end.

  - `expectedNumberOfSubjects`: Expected cumulative enrollments at trial
    end.

  - `expectedStudyDuration`: Expected study duration.

  - `criticalValues`: The input matrix of by-level critical boundaries.

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

  - `iterationNumber`, `stopStage`, `stageNumber`, `analysisTime`,
    `activeArm`, `totalAccruals`, `totalEvents`, `totalDropouts`,
    `uscore`, `vscore`, `logRankStatistic`, `reject`.

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
(sim1 = lrsim_mams(
  M = 2,
  kMax = 3,
  criticalValues = matrix(c(3.879976, 2.734557, 2.246072,
                            3.710303, 2.511427, 1.993047), 3, 2),
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  lambdas = list(log(2)/12*0.75, log(2)/12*0.75, log(2)/12),
  n = 700,
  plannedEvents = c(108, 216, 324),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 0))
#>                                                
#> Multi-arm multi-stage design for log-rank test 
#> Empirical power: 0.798                         
#> Number of active arms: 2                       
#> Number of looks: 3                             
#>                                                
#> By level critical boundaries
#>         Level 2 Level 1
#> Stage 1   3.880   3.710
#> Stage 2   2.735   2.511
#> Stage 3   2.246   1.993
#> 
#> Cumulative probability of rejection by treatment
#>         Active 1 Active 2 Overall
#> Stage 1    0.009    0.007   0.013
#> Stage 2    0.287    0.269   0.390
#> Stage 3    0.630    0.594   0.798
#> 
#> Probability of rejection by number of active arms
#>             0     1     2
#> Stage 1 0.987 0.010 0.003
#> Stage 2 0.623 0.214 0.163
#> Stage 3 0.592 0.148 0.260
#> Overall 0.202 0.372 0.426
#> 
#> Overall probability of rejection by set of active arms
#>   Set of active arms Probability of rejection
#> 1               none                    0.202
#> 2                  1                    0.204
#> 3                  2                    0.168
#> 4                1,2                    0.426
#> 
#>                            Active 1 Active 2 Control Total
#> Expected # events             128.6    128.9   151.9 409.4
#> Expected # dropouts             0.0      0.0     0.0   0.0
#> Expected # subjects           232.0    232.0   232.0 696.0
#> Expected study duration        37.9     37.9    37.9  37.9
#> 
#>                            Active 1 Active 2 Control Total
#> Number of events   Stage 1     48.2     48.5    59.8 156.5
#> Number of events   Stage 2     98.1     98.3   117.9 314.3
#> Number of events   Stage 3    151.7    151.7   172.3 475.7
#> Number of dropouts Stage 1      0.0      0.0     0.0   0.0
#> Number of dropouts Stage 2      0.0      0.0     0.0   0.0
#> Number of dropouts Stage 3      0.0      0.0     0.0   0.0
#> Number of subjects Stage 1    159.8    159.7   159.8 479.3
#> Number of subjects Stage 2    232.7    232.7   232.7 698.2
#> Number of subjects Stage 3    233.3    233.3   233.3 700.0
#> Analysis time      Stage 1     22.3     22.3    22.3  22.3
#> Analysis time      Stage 2     31.2     31.2    31.2  31.2
#> Analysis time      Stage 3     42.6     42.6    42.6  42.6
```
