# Log-Rank Test Simulation

Performs simulation for two-arm group sequential trials based on
weighted log-rank test.

## Usage

``` r
lrsim(
  kMax = 1,
  informationRates = NA,
  criticalValues = NA,
  futilityBounds = NA,
  hazardRatioH0 = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  lambda1 = NA,
  lambda2 = NA,
  gamma1 = 0,
  gamma2 = 0,
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

- informationRates:

  The information rates in terms of number of events for the
  conventional log-rank test and in terms of the actual information for
  weighted log-rank tests. Fixed prior to the trial. If left
  unspecified, it defaults to `plannedEvents / plannedEvents[kMax]` when
  `plannedEvents` is provided and to `plannedTime / plannedTime[kMax]`
  otherwise.

- criticalValues:

  Upper boundaries on the z-test statistic scale for stopping for
  efficacy.

- futilityBounds:

  Lower boundaries on the z-test statistic scale for stopping for
  futility at stages `1, ..., kMax-1`. Defaults to `rep(-8, kMax-1)` if
  left unspecified. The futility bounds are non-binding for the
  calculation of critical values.

- hazardRatioH0:

  Hazard ratio under the null hypothesis for the active treatment versus
  control. Defaults to 1 for superiority test.

- allocation1:

  Number of subjects in the active treatment group in a randomization
  block. Defaults to 1 for equal randomization.

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

- lambda1:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the active treatment group.

- lambda2:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group.

- gamma1:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the active treatment group.

- gamma2:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group.

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

  The planned cumulative total number of events at each stage.

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

An S3 class `lrsim` object with 3 components:

- `overview`: A list containing the following information:

  - `rejectPerStage`: The efficacy stopping probability by stage.

  - `futilityPerStage`: The futility stopping probability by stage.

  - `cumulativeRejection`: Cumulative efficacy stopping probability by
    stage.

  - `cumulativeFutility`: The cumulative futility stopping probability
    by stage.

  - `numberOfEvents`: The average number of events by stage.

  - `numberOfDropouts`: The average number of dropouts by stage.

  - `numberOfSubjects`: The average number of subjects by stage.

  - `analysisTime`: The average analysis time by stage.

  - `overallReject`: The overall rejection probability.

  - `expectedNumberOfEvents`: The expected number of events for the
    overall study.

  - `expectedNumberOfDropouts`: The expected number of dropouts for the
    overall study.

  - `expectedNumberOfSubjects`: The expected number of subjects for the
    overall study.

  - `expectedStudyDuration`: The expected study duration.

  - `hazardRatioH0`: Hazard ratio under the null hypothesis for the
    active treatment versus control.

  - `useEvents`: whether the analyses are planned based on the number of
    events or calendar time.

  - `numberOfIterations`: The number of simulation iterations.

  - `n`: Sample size.

  - `fixedFollowup`: Whether a fixed follow-up design is used.

  - `rho1`: The first parameter of the Fleming-Harrington family of
    weighted log-rank test. Defaults to 0 for conventional log-rank
    test.

  - `rho2`: The second parameter of the Fleming-Harrington family of
    weighted log-rank test. Defaults to 0 for conventional log-rank
    test.

  - `kMax`: The maximum number of stages.

- `sumdata`: A data frame of summary data by iteration and stage:

  - `iterationNumber`: The iteration number.

  - `eventsNotAchieved`: Whether the final target number of events is
    not achieved for the iteration.

  - `stopStage`: The stage at which the trial stops.

  - `stageNumber`: The stage number, covering all stages even if the
    trial stops at an interim look.

  - `analysisTime`: The time for the stage since trial start.

  - `accruals1`: The number of subjects enrolled at the stage for the
    treatment group.

  - `accruals2`: The number of subjects enrolled at the stage for the
    control group.

  - `totalAccruals`: The total number of subjects enrolled at the stage.

  - `events1`: The number of events at the stage for the treatment
    group.

  - `events2`: The number of events at the stage for the control group.

  - `totalEvents`: The total number of events at the stage.

  - `dropouts1`: The number of dropouts at the stage for the treatment
    group.

  - `dropouts2`: The number of dropouts at the stage for the control
    group.

  - `totalDropouts`: The total number of dropouts at the stage.

  - `uscore`: The numerator of the log-rank test statistic.

  - `vscore`: The variance of the log-rank test statistic.

  - `logRankStatistic`: The log-rank test Z-statistic.

  - `rejectPerStage`: Whether to reject the null hypothesis at the
    stage.

  - `futilityPerStage`: Whether to stop the trial for futility at the
    stage.

- `rawdata` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for selected
  replications, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stopStage`: The stage at which the trial stops.

  - `stageNumber`: The stage number, covering all stages even if the
    trial stops at an interim look.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `treatmentGroup`: The treatment group (1 or 2) for the subject.

  - `survivalTime`: The underlying survival time for the subject.

  - `dropoutTime`: The underlying dropout time for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization.

  - `event`: Whether the subject experienced the event.

  - `dropoutEvent`: Whether the subject dropped out.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Example 1: analyses based on number of events

sim1 <- lrsim(
  kMax = 2, informationRates = c(0.5, 1),
  criticalValues = c(2.797, 1.977),
  accrualIntensity = 11,
  lambda1 = 0.018, lambda2 = 0.030,
  n = 132,
  plannedEvents = c(60, 120),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)

# summary statistics
sim1
#>                                                         
#> Group-sequential design with 2 stages for log-rank test 
#> Empirical power: 0.791                                  
#> Expected # events: 106.9                                
#> Expected # dropouts: 0                                  
#> Expected # subjects: 132                                
#> Expected study duration: 94.3                           
#> n: 132, fixed follow-up: FALSE                          
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2
#> Cumulative rejection 0.2190  0.7910 
#> Cumulative futility  0.0000  0.2090 
#> Number of events     60.0    120.0  
#> Number of dropouts   0.0     0.0    
#> Number of subjects   132.0   132.0  
#> Analysis time        31.9    111.6  

# summary for each simulated data set
head(sim1$sumdata)
#>   iterationNumber stopStage eventsNotAchieved stageNumber analysisTime
#> 1               1         1             FALSE           1     33.16275
#> 2               1         1             FALSE           2    110.48736
#> 3               2         2             FALSE           1     33.22616
#> 4               2         2             FALSE           2    109.62981
#> 5               3         2             FALSE           1     28.33555
#> 6               3         2             FALSE           2    115.48439
#>   accruals1 accruals2 totalAccruals events1 events2 totalEvents dropouts1
#> 1        66        66           132      21      39          60         0
#> 2        66        66           132      56      64         120         0
#> 3        66        66           132      29      31          60         0
#> 4        66        66           132      54      66         120         0
#> 5        66        66           132      25      35          60         0
#> 6        66        66           132      58      62         120         0
#>   dropouts2 totalDropouts     uscore   vscore logRankStatistic rejectPerStage
#> 1         0             0 -11.660878 14.80310       -3.0307840           TRUE
#> 2         0             0 -18.892453 27.16956       -3.6244912           TRUE
#> 3         0             0  -1.199027 14.98331       -0.3097597          FALSE
#> 4         0             0 -12.829339 28.87298       -2.3875830           TRUE
#> 5         0             0  -6.362696 14.94729       -1.6457354          FALSE
#> 6         0             0 -11.331095 28.85624       -2.1093662           TRUE
#>   futilityPerStage
#> 1            FALSE
#> 2            FALSE
#> 3            FALSE
#> 4            FALSE
#> 5            FALSE
#> 6            FALSE

# raw data for selected replication
head(sim1$rawdata)
#>   iterationNumber stopStage stageNumber analysisTime subjectId arrivalTime
#> 1               1         1           1     33.16275         1  0.06481234
#> 2               1         1           1     33.16275         2  0.13397351
#> 3               1         1           1     33.16275         3  0.18561876
#> 4               1         1           1     33.16275         4  0.37951139
#> 5               1         1           1     33.16275         5  0.39512177
#> 6               1         1           1     33.16275         6  0.49481613
#>   stratum treatmentGroup survivalTime dropoutTime timeUnderObservation event
#> 1       1              1    23.977669         Inf            23.977669  TRUE
#> 2       1              2    46.952068         Inf            33.028772 FALSE
#> 3       1              2    21.301206         Inf            21.301206  TRUE
#> 4       1              1    32.783234         Inf            32.783234  TRUE
#> 5       1              1     6.212465         Inf             6.212465  TRUE
#> 6       1              2   122.679540         Inf            32.667930 FALSE
#>   dropoutEvent
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4        FALSE
#> 5        FALSE
#> 6        FALSE


# Example 2: analyses based on calendar time have similar power

sim2 <- lrsim(
  kMax = 2, informationRates = c(0.5, 1),
  criticalValues = c(2.797, 1.977),
  accrualIntensity = 11,
  lambda1 = 0.018, lambda2 = 0.030,
  n = 132,
  plannedTime = c(31.9, 113.2),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)

# summary statistics
sim2
#>                                                         
#> Group-sequential design with 2 stages for log-rank test 
#> Empirical power: 0.79                                   
#> Expected # events: 107.7                                
#> Expected # dropouts: 0                                  
#> Expected # subjects: 132                                
#> Expected study duration: 96.6                           
#> n: 132, fixed follow-up: FALSE                          
#> Number of simulations: 1000                             
#>                                                         
#>                      Stage 1 Stage 2
#> Cumulative rejection 0.2040  0.7900 
#> Cumulative futility  0.0000  0.2100 
#> Number of events     59.7    119.9  
#> Number of dropouts   0.0     0.0    
#> Number of subjects   132.0   132.0  
#> Analysis time        31.9    113.2  

# summary for each simulated data set
head(sim2$sumdata)
#>   iterationNumber stopStage eventsNotAchieved stageNumber analysisTime
#> 1               1         1             FALSE           1         31.9
#> 2               1         1             FALSE           2        113.2
#> 3               2         2             FALSE           1         31.9
#> 4               2         2             FALSE           2        113.2
#> 5               3         2             FALSE           1         31.9
#> 6               3         2             FALSE           2        113.2
#>   accruals1 accruals2 totalAccruals events1 events2 totalEvents dropouts1
#> 1        66        66           132      20      39          59         0
#> 2        66        66           132      57      64         121         0
#> 3        66        66           132      28      29          57         0
#> 4        66        66           132      55      66         121         0
#> 5        66        66           132      28      40          68         0
#> 6        66        66           132      57      61         118         0
#>   dropouts2 totalDropouts     uscore   vscore logRankStatistic rejectPerStage
#> 1         0             0 -11.972719 14.56264       -3.1374202           TRUE
#> 2         0             0 -18.767453 27.31019       -3.5912283           TRUE
#> 3         0             0  -0.632938 14.22525       -0.1678152          FALSE
#> 4         0             0 -12.970248 28.77936       -2.4177294           TRUE
#> 5         0             0  -8.140083 16.88593       -1.9809173          FALSE
#> 6         0             0 -10.971205 28.42533       -2.0577926           TRUE
#>   futilityPerStage
#> 1            FALSE
#> 2            FALSE
#> 3            FALSE
#> 4            FALSE
#> 5            FALSE
#> 6            FALSE
```
