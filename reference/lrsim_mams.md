# Log-Rank Test Simulation for Multi-Arm Multi-Stage Design

Simulate a multi-arm multi-stage design using a weighted log-rank test.
Analyses can be triggered either by the cumulative number of events
(combined for an active arm and the common control) or by pre-specified
calendar times.

## Usage

``` r
lrsim_mams(
  M = 2,
  kMax = 1,
  criticalValues = NULL,
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

- futilityBounds:

  Numeric vector of length `kMax - 1` giving the futility boundaries on
  the max-Z scale for the first `kMax - 1` analyses. At an interim look,
  the study stops for futility if all active treatment arms cross the
  futility boundary. At the final look, the study is counted as stopping
  for futility if none of the active treatment arms can be rejected. If
  omitted, no interim futility stopping is applied.

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

  - `overallReject`: Overall probability of rejecting the null by trial
    end.

  - `overallFutility`: Overall probability of stopping for futility by
    trial end.

  - `rejectPerStage`: Probability of rejecting the null for each active
    arm at each stage.

  - `futilityPerStage`: Probability of futility stopping for each arm at
    each stage.

  - `cumulativeRejection`: Cumulative probability of rejection for each
    active arm by stage.

  - `cumulativeFutility`: Cumulative probability of futility stopping
    for each active arm by stage.

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

  - `futilityBounds`: The input vector of futility boundaries.

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
    `uscore`, `vscore`, `logRankStatistic`, `reject`, `futility`.

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
(sim1 <- lrsim_mams(
  M = 2,
  kMax = 3,
  criticalValues = matrix(c(3.880, 2.747, 2.275,
                            3.710, 2.511, 1.993), 3, 2),
  futilityBounds = c(0.074, 1.207),
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  lambdas = list(log(2)/12*0.5, log(2)/12*0.75, log(2)/12),
  n = 700,
  plannedEvents = c(36, 72, 108),
  maxNumberOfIterations = 10000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 0))
#>                                                
#> Multi-arm multi-stage design for log-rank test 
#> Empirical power: 0.8963                        
#> Number of active arms: 2                       
#> Number of looks: 3                             
#>                                                
#> By level critical boundaries
#>         Level 2 Level 1
#> Stage 1   3.880   3.710
#> Stage 2   2.747   2.511
#> Stage 3   2.275   1.993
#> 
#> Cumulative probability of rejection or futility by treatment
#>         Reject Active 1 Reject Active 2 Overall Rejection Futility
#> Stage 1          0.0244          0.0011            0.0248   0.0139
#> Stage 2          0.5596          0.0928            0.5670   0.0422
#> Stage 3          0.8870          0.1802            0.8963   0.1037
#> 
#> Detailed probability of trial termination at each look
#>         Reject 1 Active Reject 2 Actives Overall Rejection Futility Continue
#> Stage 1          0.0241           0.0007            0.0248   0.0139   0.9613
#> Stage 2          0.4575           0.0847            0.5422   0.0283   0.3908
#> Stage 3          0.2438           0.0855            0.3293   0.0615   0.0000
#> Total            0.7254           0.1709            0.8963   0.1037       NA
#> 
#> Overall probability of rejection by set of active arms
#>   Set of active arms Probability of rejection
#> 1               none                   0.1037
#> 2                  1                   0.7161
#> 3                  2                   0.0093
#> 4                1,2                   0.1709
#> 
#>                            Active 1 Active 2 Control Total
#> Expected # events              31.0     43.1    53.7 127.8
#> Expected # dropouts             0.0      0.0     0.0   0.0
#> Expected # subjects           148.4    148.4   148.4 445.2
#> Expected study duration        21.0     21.0    21.0  21.0
#> 
#>                            Active 1 Active 2 Control Total
#> Number of events   Stage 1     12.9     18.2    23.1  54.2
#> Number of events   Stage 2     26.3     36.7    45.7 108.7
#> Number of events   Stage 3     42.2     55.0    65.8 163.0
#> Number of dropouts Stage 1      0.0      0.0     0.0   0.0
#> Number of dropouts Stage 2      0.0      0.0     0.0   0.0
#> Number of dropouts Stage 3      0.0      0.0     0.0   0.0
#> Number of subjects Stage 1     89.4     89.4    89.4 268.3
#> Number of subjects Stage 2    135.8    135.8   135.8 407.5
#> Number of subjects Stage 3    172.0    172.0   172.0 516.1
#> Analysis time      Stage 1     14.7     14.7    14.7  14.7
#> Analysis time      Stage 2     19.7     19.7    19.7  19.7
#> Analysis time      Stage 3     23.6     23.6    23.6  23.6
```
