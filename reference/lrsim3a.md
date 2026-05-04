# Log-Rank Test Simulation for Three Arms

Performs simulation for three-arm group sequential trials based on
weighted log-rank test. The looks are driven by the total number of
events in Arm A and Arm C combined. Alternatively, the analyses can be
planned to occur at specified calendar times.

## Usage

``` r
lrsim3a(
  kMax = 1,
  hazardRatioH013 = 1,
  hazardRatioH023 = 1,
  hazardRatioH012 = 1,
  allocation1 = 1,
  allocation2 = 1,
  allocation3 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  lambda1 = NA,
  lambda2 = NA,
  lambda3 = NA,
  gamma1 = 0,
  gamma2 = 0,
  gamma3 = 0,
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

- hazardRatioH013:

  Hazard ratio under the null hypothesis for arm 1 versus arm 3.
  Defaults to 1 for superiority test.

- hazardRatioH023:

  Hazard ratio under the null hypothesis for arm 2 versus arm 3.
  Defaults to 1 for superiority test.

- hazardRatioH012:

  Hazard ratio under the null hypothesis for arm 1 versus arm 2.
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

- lambda1:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 1.

- lambda2:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 2.

- lambda3:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for arm 3.

- gamma1:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 1.

- gamma2:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 2.

- gamma3:

  The hazard rate for exponential dropout. A vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for arm 3.

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

  The planned cumulative total number of events at Look 1 to Look `kMax`
  for Arms A and C combined.

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

  - `uscore13`: The log-rank test score statistic comparing the active
    treatment 1 to the control.

  - `vscore13`: The log-rank test variance statistic comparing the
    active treatment 1 to the control.

  - `logRankStatistic13`: The log-rank test Z-statistic comparing the
    active treatment 1 to the control.

  - `uscore23`: The log-rank test score statistic comparing the active
    treatment 2 to the control.

  - `vscore23`: The log-rank test variance statistic comparing the
    active treatment 2 to the control.

  - `logRankStatistic23`: The log-rank test Z-statistic comparing the
    active treatment 2 to the control.

  - `uscore12`: The log-rank test score statistic comparing the active
    treatment 1 to the active treatment 2.

  - `vscore12`: The log-rank test variance statistic comparing the
    active treatment 1 to the active treatment 2.

  - `logRankStatistic12`: The log-rank test Z-statistic comparing the
    active treatment 1 to the active treatment 2.

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

  - `survivalTime`: The underlying survival time for the subject.

  - `dropoutTime`: The underlying dropout time for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization for the subject.

  - `event`: Whether the subject experienced the event.

  - `dropoutEvent`: Whether the subject dropped out.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
sim1 <- lrsim3a(
  kMax = 3,
  allocation1 = 2,
  allocation2 = 2,
  allocation3 = 1,
  accrualTime = c(0, 8),
  accrualIntensity = c(10, 28),
  piecewiseSurvivalTime = 0,
  lambda1 = log(2)/12*0.60,
  lambda2 = log(2)/12*0.70,
  lambda3 = log(2)/12,
  n = 700,
  plannedEvents = c(186, 259, 295),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)

head(sim1$sumdata)
#>   iterationNumber eventsNotAchieved stageNumber analysisTime accruals1
#> 1               1             FALSE           1     32.50989       280
#> 2               1             FALSE           2     41.81493       280
#> 3               1             FALSE           3     48.34845       280
#> 4               2             FALSE           1     34.76899       280
#> 5               2             FALSE           2     44.51001       280
#> 6               2             FALSE           3     51.77071       280
#>   accruals2 accruals3 totalAccruals events1 events2 events3 totalEvents
#> 1       280       140           700     112     118      74         304
#> 2       280       140           700     163     169      96         428
#> 3       280       140           700     182     192     113         487
#> 4       280       140           700     104     136      82         322
#> 5       280       140           700     153     182     106         441
#> 6       280       140           700     174     207     121         502
#>   dropouts1 dropouts2 dropouts3 totDropouts  uscore13 vscore13
#> 1         0         0         0           0 -16.94859 39.39278
#> 2         0         0         0           0 -18.59502 53.93141
#> 3         0         0         0           0 -26.47463 60.69003
#> 4         0         0         0           0 -28.22168 37.95361
#> 5         0         0         0           0 -35.27386 51.00417
#> 6         0         0         0           0 -44.13876 56.08485
#>   logRankStatistic13  uscore23 vscore23 logRankStatistic23   uscore12 vscore12
#> 1          -2.700383 -15.39159 40.57090          -2.416441  -2.896980 57.44915
#> 2          -2.532070 -15.58195 55.62798          -2.089177  -5.322203 82.51706
#> 3          -3.398374 -22.62257 63.07148          -2.848560  -7.654443 93.31583
#> 4          -4.580957 -12.02085 47.45829          -1.744935 -23.881566 59.60921
#> 5          -4.939127 -16.52496 61.55452          -2.106252 -28.451270 82.94223
#> 6          -5.893826 -21.85525 68.78296          -2.635210 -34.578948 94.04171
#>   logRankStatistic12
#> 1         -0.3822115
#> 2         -0.5858947
#> 3         -0.7923844
#> 4         -3.0931864
#> 5         -3.1240199
#> 6         -3.5657552
head(sim1$rawdata)
#>   iterationNumber stageNumber analysisTime subjectId arrivalTime stratum
#> 1               1           1     32.50989         1  0.07129358       1
#> 2               1           1     32.50989         2  0.14737086       1
#> 3               1           1     32.50989         3  0.20418064       1
#> 4               1           1     32.50989         4  0.41746253       1
#> 5               1           1     32.50989         5  0.43463395       1
#> 6               1           1     32.50989         6  0.54429775       1
#>   treatmentGroup survivalTime dropoutTime timeUnderObservation event
#> 1              1    12.453287         Inf            12.453287  TRUE
#> 2              1    40.642510         Inf            32.362519 FALSE
#> 3              2    15.804589         Inf            15.804589  TRUE
#> 4              3    10.215981         Inf            10.215981  TRUE
#> 5              2     2.765631         Inf             2.765631  TRUE
#> 6              3    63.716099         Inf            31.965593 FALSE
#>   dropoutEvent
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4        FALSE
#> 5        FALSE
#> 6        FALSE
```
