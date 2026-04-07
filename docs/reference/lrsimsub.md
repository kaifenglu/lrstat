# Log-Rank Test Simulation for Enrichment Design

Performs simulation for two-arm group sequential trials based on
weighted log-rank test for a biomarker enrichment design. The looks are
either driven by the total number of events in the ITT population or the
biomarker positive sub population. Alternatively, the analyses can be
planned to occur at specified calendar times.

## Usage

``` r
lrsimsub(
  kMax = 1,
  kMaxitt = 1,
  hazardRatioH0itt = 1,
  hazardRatioH0pos = 1,
  hazardRatioH0neg = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  p_pos = NA,
  lambda1itt = NA,
  lambda2itt = NA,
  lambda1pos = NA,
  lambda2pos = NA,
  gamma1itt = 0,
  gamma2itt = 0,
  gamma1pos = 0,
  gamma2pos = 0,
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

- kMaxitt:

  Number of stages with timing determined by events in the ITT
  population. Ranges from 0 (none) to `kMax`.

- hazardRatioH0itt:

  Hazard ratio under the null hypothesis for the ITT population.
  Defaults to 1 for superiority test.

- hazardRatioH0pos:

  Hazard ratio under the null hypothesis for the biomarker positive sub
  population. Defaults to 1 for superiority test.

- hazardRatioH0neg:

  Hazard ratio under the null hypothesis for the biomarker negative sub
  population. Defaults to 1 for superiority test.

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

- p_pos:

  The prevalence of the biomarker positive sub population in each
  stratum.

- lambda1itt:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the treatment group in the ITT population.

- lambda2itt:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group in the ITT population.

- lambda1pos:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the treatment group in the biomarker positive sub
  population.

- lambda2pos:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group in the biomarker positive sub
  population.

- gamma1itt:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the treatment group in the ITT population.

- gamma2itt:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group in the ITT population.

- gamma1pos:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the treatment group in the biomarker positive sub population.

- gamma2pos:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group in the biomarker positive sub population.

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

  The planned cumulative total number events in the ITT population at
  Look 1 to Look `kMaxitt` and the planned cumulative total number of
  events at Look `kMaxitt+1` to Look `kMax` in the biomarker positive
  sub population.

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

  - `population`: The population ("ITT", "Biomarker Positive",
    "Biomarker Negative") under consideration.

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

  - `logRankStatistic`: The log-rank test Z-statistic for the
    population.

- `rawdata` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for selected
  replications, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage under consideration.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `biomarker`: The biomarker status for the subject (1 for positive, 0
    for negative).

  - `treatmentGroup`: The treatment group (1 or 2) for the subject.

  - `survivalTime`: The underlying survival time for the subject.

  - `dropoutTime`: The underlying dropout time for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization for the subject.

  - `event`: Whether the subject experienced an event.

  - `dropoutEvent`: Whether the subject dropped out.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
sim1 = lrsimsub(
  kMax = 2,
  kMaxitt = 2,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = seq(0,9),
  accrualIntensity = c(seq(10,70,10),rep(70,3)),
  piecewiseSurvivalTime = c(0,12,24),
  p_pos = 0.6,
  lambda1itt = c(0.00256, 0.00383, 0.00700),
  lambda2itt = c(0.00427, 0.00638, 0.01167),
  lambda1pos = c(0.00299, 0.00430, 0.01064),
  lambda2pos = c(0.00516, 0.00741, 0.01835),
  gamma1itt = -log(1-0.04)/12,
  gamma2itt = -log(1-0.04)/12,
  gamma1pos = -log(1-0.04)/12,
  gamma2pos = -log(1-0.04)/12,
  n = 500,
  plannedEvents = c(108,144),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)

head(sim1$sumdata)
#>   iterNumber events1NotAchieved events2NotAchieved stageNumber analysisTime
#> 1          1              FALSE              FALSE           1     47.20050
#> 2          1              FALSE              FALSE           1     47.20050
#> 3          1              FALSE              FALSE           1     47.20050
#> 4          1              FALSE              FALSE           2     59.14058
#> 5          1              FALSE              FALSE           2     59.14058
#> 6          1              FALSE              FALSE           2     59.14058
#>   population accruals1 accruals2 totalAccruals events1 events2 totalEvents
#> 1        ITT       250       250           500      38      70         108
#> 2 Biomarker+       143       136           279      32      45          77
#> 3 Biomarker-       107       114           221       6      25          31
#> 4        ITT       250       250           500      56      88         144
#> 5 Biomarker+       143       136           279      49      62         111
#> 6 Biomarker-       107       114           221       7      26          33
#>   dropouts1 dropouts2 totalDropouts     uscore    vscore logRankStatistic
#> 1        29        39            68 -19.822223 26.812336        -3.828115
#> 2        16        23            39 -11.136324 18.868887        -2.563709
#> 3        13        16            29  -9.566429  7.687120        -3.450389
#> 4        36        45            81 -22.525041 35.622174        -3.774030
#> 5        20        28            48 -14.240702 27.082819        -2.736431
#> 6        16        17            33 -10.008313  8.228709        -3.488952
head(sim1$rawdata)
#>   iterationNumber stageNumber analysisTime subjectId arrivalTime stratum
#> 1               1           1      47.2005         1  0.07129358       1
#> 2               1           1      47.2005         2  0.15135848       1
#> 3               1           1      47.2005         3  0.25417963       1
#> 4               1           1      47.2005         4  0.31318945       1
#> 5               1           1      47.2005         5  0.48171673       1
#> 6               1           1      47.2005         6  0.74681987       1
#>   biomarker treatmentGroup survivalTime dropoutTime timeUnderObservation event
#> 1      TRUE              1     14.73791   223.63615             14.73791  TRUE
#> 2     FALSE              2    140.84362   269.95567             47.04914 FALSE
#> 3      TRUE              2    154.76417   802.30215             46.94632 FALSE
#> 4      TRUE              1     32.82033    32.87180             32.82033  TRUE
#> 5     FALSE              2    521.59369   235.55130             46.71878 FALSE
#> 6      TRUE              1     53.80430    31.18759             31.18759 FALSE
#>   dropoutEvent
#> 1        FALSE
#> 2        FALSE
#> 3        FALSE
#> 4        FALSE
#> 5        FALSE
#> 6         TRUE
```
