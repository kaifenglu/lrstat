# Log-Rank Test Simulation for Two Endpoints (PFS and OS) and Three Arms

Performs simulation for two-endpoint (PFS and OS) three-arm group
sequential trials based on weighted log-rank test. The first `kMaxpfs`
looks are driven by the total number of PFS events in Arm A and Arm C
combined, and the subsequent looks are driven by the total number of OS
events in Arm A and Arm C combined. Alternatively, the analyses can be
planned to occur at specified calendar times.

## Usage

``` r
lrsim2e3a(
  kMax = 1,
  kMaxpfs = 1,
  hazardRatioH013pfs = 1,
  hazardRatioH023pfs = 1,
  hazardRatioH012pfs = 1,
  hazardRatioH013os = 1,
  hazardRatioH023os = 1,
  hazardRatioH012os = 1,
  allocation1 = 1,
  allocation2 = 1,
  allocation3 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  rho_pd_os = 0,
  lambda1pfs = NA,
  lambda2pfs = NA,
  lambda3pfs = NA,
  lambda1os = NA,
  lambda2os = NA,
  lambda3os = NA,
  gamma1pfs = 0,
  gamma2pfs = 0,
  gamma3pfs = 0,
  gamma1os = 0,
  gamma2os = 0,
  gamma3os = 0,
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

- kMax:

  The maximum number of stages.

- kMaxpfs:

  Number of stages with timing determined by PFS events. Ranges from 0
  (none) to `kMax`.

- hazardRatioH013pfs:

  Hazard ratio under the null hypothesis for arm 1 vs arm 3 for PFS.
  Defaults to 1 for superiority test.

- hazardRatioH023pfs:

  Hazard ratio under the null hypothesis for arm 2 vs arm 3 for PFS.
  Defaults to 1 for superiority test.

- hazardRatioH012pfs:

  Hazard ratio under the null hypothesis for arm 1 vs arm 2 for PFS.
  Defaults to 1 for superiority test.

- hazardRatioH013os:

  Hazard ratio under the null hypothesis for arm 1 vs arm 3 for OS.
  Defaults to 1 for superiority test.

- hazardRatioH023os:

  Hazard ratio under the null hypothesis for arm 2 vs arm 3 for OS.
  Defaults to 1 for superiority test.

- hazardRatioH012os:

  Hazard ratio under the null hypothesis for arm 1 vs arm 2 for OS.
  Defaults to 1 for superiority test.

- allocation1:

  Number of subjects in Arm A in a randomization block. Defaults to 1
  for equal randomization.

- allocation2:

  Number of subjects in Arm B in a randomization block. Defaults to 1
  for equal randomization.

- allocation3:

  Number of subjects in Arm C in a randomization block. Defaults to 1
  for equal randomization.

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

- rho_pd_os:

  The correlation coefficient for the standard bivariate normal random
  variables used to generate time to disease progression and time to
  death using the inverse CDF method.

- lambda1pfs:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 1 and PFS.

- lambda2pfs:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 2 and PFS.

- lambda3pfs:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 3 and PFS.

- lambda1os:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 1 and OS.

- lambda2os:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 2 and OS.

- lambda3os:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 3 and OS.

- gamma1pfs:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 1 and PFS.

- gamma2pfs:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 2 and PFS.

- gamma3pfs:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 3 and PFS.

- gamma1os:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 1 and OS.

- gamma2os:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 2 and OS.

- gamma3os:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 3 and OS.

- n:

  Sample size.

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

  The planned cumulative total number of PFS events at Look 1 to Look
  `kMaxpfs` for Arms A and C combined and the planned cumulative total
  number of OS events at Look `kMaxpfs+1` to Look `kMax` for Arms A and
  C combined.

- plannedTime:

  The calendar times for the analyses. To use calendar time to plan the
  analyses, `plannedEvents` should be missing.

- maxNumberOfIterations:

  The number of simulation iterations. Defaults to 1000.

- maxNumberOfRawDatasetsPerStage:

  The number of raw datasets per stage to extract.

- seed:

  The seed to reproduce the simulation results.

- nthreads:

  The number of threads to use in simulations (0 means the default
  RcppParallel behavior).

## Value

A list with 2 components:

- `sumdata`: A data frame of summary data by iteration and stage:

  - `iterationNumber`: The iteration number.

  - `eventsNotAchieved`: Whether the target number of events is not
    achieved for the iteration.

  - `stageNumber`: The stage number, covering all stages even if the
    trial stops at an interim look.

  - `analysisTime`: The time for the stage since trial start.

  - `accruals1`: The number of subjects enrolled at the stage for the
    active treatment 1 group.

  - `accruals2`: The number of subjects enrolled at the stage for the
    active treatment 2 group.

  - `accruals3`: The number of subjects enrolled at the stage for the
    control group.

  - `totalAccruals`: The total number of subjects enrolled at the stage.

  - `endpoint`: The endpoint (1 for PFS or 2 for OS) under
    consideration.

  - `events1`: The number of events at the stage for the active
    treatment 1 group.

  - `events2`: The number of events at the stage for the active
    treatment 2 group.

  - `events3`: The number of events at the stage for the control group.

  - `totalEvents`: The total number of events at the stage.

  - `dropouts1`: The number of dropouts at the stage for the active
    treatment 1 group.

  - `dropouts2`: The number of dropouts at the stage for the active
    treatment 2 group.

  - `dropouts3`: The number of dropouts at the stage for the control
    group.

  - `totalDropouts`: The total number of dropouts at the stage.

  - `logRankStatistic13`: The log-rank test Z-statistic comparing the
    active treatment 1 to the control for the endpoint.

  - `logRankStatistic23`: The log-rank test Z-statistic comparing the
    active treatment 2 to the control for the endpoint.

  - `logRankStatistic12`: The log-rank test Z-statistic comparing the
    active treatment 1 to the active treatment 2 for the endpoint.

- `rawdata` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for selected
  replications, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage under consideration.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `treatmentGroup`: The treatment group (1, 2, or 3) for the subject.

  - `endpoint`: The endpoint (1 for PFS or 2 for OS) under
    consideration.

  - `survivalTime`: The underlying survival time for the event endpoint
    for the subject.

  - `dropoutTime`: The underlying dropout time for the event endpoint
    for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization for the event endpoint for the subject.

  - `event`: Whether the subject experienced the event endpoint.

  - `dropoutEvent`: Whether the subject dropped out for the endpoint.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
sim1 <- lrsim2e3a(
  kMax = 3,
  kMaxpfs = 2,
  allocation1 = 2,
  allocation2 = 2,
  allocation3 = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  rho_pd_os = 0,
  lambda1pfs = log(2)/12*0.60,
  lambda2pfs = log(2)/12*0.70,
  lambda3pfs = log(2)/12,
  lambda1os = log(2)/30*0.65,
  lambda2os = log(2)/30*0.75,
  lambda3os = log(2)/30,
  n = 700,
  plannedEvents = c(186, 259, 183),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)

head(sim1$sumdata)
#>   iterationNumber events1NotAchieved events2NotAchieved stageNumber
#> 1               1              FALSE              FALSE           1
#> 2               1              FALSE              FALSE           1
#> 3               1              FALSE              FALSE           2
#> 4               1              FALSE              FALSE           2
#> 5               1              FALSE              FALSE           3
#> 6               1              FALSE              FALSE           3
#>   analysisTime accruals1 accruals2 accruals3 totalAccruals endpoint events1
#> 1     32.85250       280       280       140           700      PFS     111
#> 2     32.85250       280       280       140           700       OS      65
#> 3     42.89920       280       280       140           700      PFS     159
#> 4     42.89920       280       280       140           700       OS      96
#> 5     49.48199       280       280       140           700      PFS     177
#> 6     49.48199       280       280       140           700       OS     110
#>   events2 events3 totalEvents dropouts1 dropouts2 dropouts3 totalDropouts
#> 1     125      75         311         0         0         0             0
#> 2      69      40         174         0         0         0             0
#> 3     187     100         446         0         0         0             0
#> 4     111      63         270         0         0         0             0
#> 5     209     118         504         0         0         0             0
#> 6     132      73         315         0         0         0             0
#>     uscore13 vscore13 logRankStatistic13   uscore23 vscore23 logRankStatistic23
#> 1 -19.813095 38.66058          -3.186531 -12.844410 42.65135          -1.966744
#> 2  -6.266267 22.87605          -1.310142  -5.015437 23.66224          -1.031053
#> 3 -25.924031 52.54320          -3.576385 -12.152275 60.83790          -1.558011
#> 4 -12.775991 34.30227          -2.181389  -7.099940 37.91929          -1.152987
#> 5 -35.700291 58.91287          -4.651218 -18.753487 68.82455          -2.260530
#> 6 -15.938713 39.18408          -2.546235  -7.314897 44.61121          -1.095182
#>     uscore12 vscore12 logRankStatistic12
#> 1 -11.937995 58.62876          -1.559107
#> 2  -2.625340 33.49406          -0.453630
#> 3 -24.461124 85.58818          -2.644049
#> 4  -9.406536 51.68518          -1.308419
#> 5 -28.803833 95.43596          -2.948453
#> 6 -13.890340 60.39008          -1.787434
head(sim1$rawdata)
#>   iterationNumber stageNumber analysisTime subjectId arrivalTime stratum
#> 1               1           1      32.8525         1  0.07129358       1
#> 2               1           1      32.8525         1  0.07129358       1
#> 3               1           1      32.8525         2  0.08504918       1
#> 4               1           1      32.8525         2  0.08504918       1
#> 5               1           1      32.8525         3  0.11887660       1
#> 6               1           1      32.8525         3  0.11887660       1
#>   treatmentGroup endpoint survivalTime dropoutTime timeUnderObservation event
#> 1              1      PFS    51.581804         Inf            32.781207 FALSE
#> 2              1       OS    52.812341         Inf            32.781207 FALSE
#> 3              1      PFS    81.221695         Inf            32.767452 FALSE
#> 4              1       OS   136.515641         Inf            32.767452 FALSE
#> 5              3      PFS     3.331084         Inf             3.331084  TRUE
#> 6              3       OS    18.301343         Inf            18.301343  TRUE
#>   dropoutEvent
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4        FALSE
#> 5        FALSE
#> 6        FALSE
```
