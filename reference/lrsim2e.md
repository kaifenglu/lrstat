# Log-Rank Test Simulation for PFS and OS Endpoints

Performs simulation for two-endpoint (PFS and OS) two-arm group
sequential trials based on weighted log-rank test. The first `kMaxpfs`
looks are driven by the total number of PFS events in two arms combined,
and the subsequent looks are driven by the total number of OS events in
two arms combined. Alternatively, the analyses can be planned to occur
at specified calendar times.

## Usage

``` r
lrsim2e(
  kMax = 1,
  kMaxpfs = 1,
  hazardRatioH0pfs = 1,
  hazardRatioH0os = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  rho_pd_os = 0,
  lambda1pfs = NA,
  lambda2pfs = NA,
  lambda1os = NA,
  lambda2os = NA,
  gamma1pfs = 0,
  gamma2pfs = 0,
  gamma1os = 0,
  gamma2os = 0,
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

- hazardRatioH0pfs:

  Hazard ratio under the null hypothesis for the active treatment vs
  control for PFS. Defaults to 1 for superiority test.

- hazardRatioH0os:

  Hazard ratio under the null hypothesis for the active treatment vs
  control for OS. Defaults to 1 for superiority test.

- allocation1:

  Number of subjects in the treatment group in a randomization block.
  Defaults to 1 for equal randomization.

- allocation2:

  Number of subjects in the control group in a randomization block.
  Defaults to 1 for equal randomization.

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
  variables used to generate time to disease progression (PD) and time
  to death using the inverse CDF method.

- lambda1pfs:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the treatment group and PFS.

- lambda2pfs:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group and PFS.

- lambda1os:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the treatment group and OS.

- lambda2os:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group and OS.

- gamma1pfs:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the treatment group and PFS.

- gamma2pfs:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group and PFS.

- gamma1os:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the treatment group and OS.

- gamma2os:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group and OS.

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
  `kMaxpfs` and the planned cumulative total number of OS events at Look
  `kMaxpfs+1` to Look `kMax`.

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
    treatment group.

  - `accruals2`: The number of subjects enrolled at the stage for the
    control group.

  - `totalAccruals`: The total number of subjects enrolled at the stage.

  - `endpoint`: The endpoint (1 for PFS or 2 for OS) under
    consideration.

  - `events1`: The number of events at the stage for the treatment
    group.

  - `events2`: The number of events at the stage for the control group.

  - `totalEvents`: The total number of events at the stage.

  - `dropouts1`: The number of dropouts at the stage for the treatment
    group.

  - `dropouts2`: The number of dropouts at the stage for the control
    group.

  - `totalDropouts`: The total number of dropouts at the stage.

  - `uscore`: The numerator of the log-rank test statistic for the
    endpoint.

  - `vscore`: The variance of the log-rank test statistic for the
    endpoint.

  - `logRankStatistic`: The log-rank test Z-statistic for the endpoint.

- `rawdata` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for selected
  replications, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage under consideration.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `treatmentGroup`: The treatment group (1 or 2) for the subject.

  - `endpoint`: The endpoint (1 for PFS or 2 for OS) under consideration
    for the row. Each subject will have two rows, one for each endpoint.

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
sim1 <- lrsim2e(
  kMax = 3,
  kMaxpfs = 2,
  allocation1 = 2,
  allocation2 = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  rho_pd_os = 0,
  lambda1pfs = log(2)/12*0.60,
  lambda2pfs = log(2)/12,
  lambda1os = log(2)/30*0.65,
  lambda2os = log(2)/30,
  n = 420,
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
#>   analysisTime accruals1 accruals2 totalAccruals endpoint events1 events2
#> 1     27.70630       280       140           420      PFS     103      83
#> 2     27.70630       280       140           420       OS      51      49
#> 3     37.52269       280       140           420      PFS     151     108
#> 4     37.52269       280       140           420       OS      79      65
#> 5     47.22055       280       140           420      PFS     186     120
#> 6     47.22055       280       140           420       OS     106      77
#>   totalEvents dropouts1 dropouts2 totalDropouts    uscore   vscore
#> 1         186         0         0             0 -28.79721 38.19448
#> 2         100         0         0             0 -17.97413 21.35703
#> 3         259         0         0             0 -37.96130 50.38427
#> 4         144         0         0             0 -22.15222 29.94302
#> 5         306         0         0             0 -41.09552 57.62563
#> 6         183         0         0             0 -23.62725 37.66037
#>   logRankStatistic
#> 1        -4.659615
#> 2        -3.889354
#> 3        -5.348027
#> 4        -4.048271
#> 5        -5.413606
#> 6        -3.850090
head(sim1$rawdata)
#>   iterationNumber stageNumber analysisTime subjectId arrivalTime stratum
#> 1               1           1      27.7063         1  0.07129358       1
#> 2               1           1      27.7063         1  0.07129358       1
#> 3               1           1      27.7063         2  0.08504918       1
#> 4               1           1      27.7063         2  0.08504918       1
#> 5               1           1      27.7063         3  0.11887660       1
#> 6               1           1      27.7063         3  0.11887660       1
#>   treatmentGroup endpoint survivalTime dropoutTime timeUnderObservation event
#> 1              1      PFS    51.581804         Inf            27.635011 FALSE
#> 2              1       OS    52.812341         Inf            27.635011 FALSE
#> 3              1      PFS    81.221695         Inf            27.621255 FALSE
#> 4              1       OS   136.515641         Inf            27.621255 FALSE
#> 5              2      PFS     3.331084         Inf             3.331084  TRUE
#> 6              2       OS    18.301343         Inf            18.301343  TRUE
#>   dropoutEvent
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4        FALSE
#> 5        FALSE
#> 6        FALSE
```
