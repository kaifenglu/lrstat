# Simulation for a Binary and a Time-to-Event Endpoint in Group Sequential Trials

Performs simulation for two-endpoint two-arm group sequential trials.

- Endpoint 1: Binary endpoint, analyzed using the Mantel-Haenszel test
  for risk difference.

- Endpoint 2: Time-to-event endpoint, analyzed using the log-rank test
  for treatment effect.

The analysis times for the binary endpoint are based on calendar times,
while the time-to-event analyses are triggered by reaching the
pre-specified number of events. The binary endpoint is assessed at the
first post-treatment follow-up visit (PTFU1).

## Usage

``` r
binary_tte_sim(
  kMax1 = 1,
  kMax2 = 1,
  riskDiffH0 = 0,
  hazardRatioH0 = 1,
  allocation1 = 1,
  allocation2 = 1,
  accrualTime = 0,
  accrualIntensity = NA,
  piecewiseSurvivalTime = 0,
  stratumFraction = 1,
  globalOddsRatio = NA,
  pi1 = NA,
  pi2 = NA,
  lambda1 = NA,
  lambda2 = NA,
  gamma1 = NA,
  gamma2 = NA,
  delta1 = NA,
  delta2 = NA,
  upper1 = NA,
  upper2 = NA,
  n = NA,
  plannedTime = NA,
  plannedEvents = NA,
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 0,
  seed = 0,
  nthreads = 0
)
```

## Arguments

- kMax1:

  Number of stages for the binary endpoint.

- kMax2:

  Number of stages for the time-to-event endpoint.

- riskDiffH0:

  Risk difference under the null hypothesis for the binary endpoint.

- hazardRatioH0:

  Hazard ratio under the null hypothesis for the time-to-event endpoint.

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

- globalOddsRatio:

  Global odds ratio of the Plackett copula linking the two endpoints.

- pi1:

  Response probabilities by stratum for the treatment group for the
  binary endpoint.

- pi2:

  Response probabilities by stratum for the control group for the binary
  endpoint.

- lambda1:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the treatment group for the time-to-event endpoint.

- lambda2:

  A vector of hazard rates for the event in each analysis time interval
  by stratum for the control group for the time-to-event endpoint.

- gamma1:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the treatment group.

- gamma2:

  The hazard rate for exponential dropout, a vector of hazard rates for
  piecewise exponential dropout applicable for all strata, or a vector
  of hazard rates for dropout in each analysis time interval by stratum
  for the control group.

- delta1:

  The hazard rate for exponential treatment discontinuation, a vector of
  hazard rates for piecewise exponential treatment discontinuation
  applicable for all strata, or a vector of hazard rates for treatment
  discontinuation in each analysis time interval by stratum for the
  treatment group for the binary endpoint.

- delta2:

  The hazard rate for exponential treatment discontinuation, a vector of
  hazard rates for piecewise exponential treatment discontinuation
  applicable for all strata, or a vector of hazard rates for treatment
  discontinuation in each analysis time interval by stratum for the
  control group for the binary endpoint.

- upper1:

  Maximim protocol-specified treatment duration for the treatment group.

- upper2:

  Maximum protocol-specified treatment duration for the control group.

- n:

  Sample size.

- plannedTime:

  Calendar times for the analyses of the binary endpoint.

- plannedEvents:

  Target cumulative number of events for the time-to-event analyses.

- maxNumberOfIterations:

  Number of simulation iterations to perform.

- maxNumberOfRawDatasetsPerStage:

  Number of subject-level datasets to retain per stage. Set to 0 to skip
  raw data saving.

- seed:

  The seed to reproduce the simulation results.

- nthreads:

  The number of threads to use in simulations (0 means the default
  RcppParallel behavior).

## Value

A list with 4 components:

- `sumdataBIN`: A data frame of summary data by iteration and stage for
  the binary endpoint:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage number, covering all stages even if the
    trial stops at an interim look.

  - `analysisTime`: The time for the stage since trial start.

  - `accruals1`: The number of subjects enrolled at the stage for the
    treatment group.

  - `accruals2`: The number of subjects enrolled at the stage for the
    control group.

  - `totalAccruals`: The total number of subjects enrolled at the stage.

  - `source1`: The total number of subjects with response status
    determined by the underlying latent response variable.

  - `source2`: The total number of subjects with response status
    (non-responder) determined by experiencing the event for the
    time-to-event endpoint.

  - `source3`: The total number of subjects with response status
    (non-responder) determined by dropping out prior to the PTFU1 visit.

  - `n1`: The number of subjects included in the analysis of the binary
    endpoint for the treatment group.

  - `n2`: The number of subjects included in the analysis of the binary
    endpoint for the control group.

  - `n`: The total number of subjects included in the analysis of the
    binary endpoint at the stage.

  - `y1`: The number of responders for the binary endpoint in the
    treatment group.

  - `y2`: The number of responders for the binary endpoint in the
    control group.

  - `y`: The total number of responders for the binary endpoint at the
    stage.

  - `riskDiff`: The estimated risk difference for the binary endpoint.

  - `seRiskDiff`: The standard error for risk difference based on the
    Sato approximation.

  - `mhStatistic`: The Mantel-Haenszel test Z-statistic for the binary
    endpoint.

- `sumdataTTE`: A data frame of summary data by iteration and stage for
  the time-to-event endpoint:

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
    time-to-event endpoint.

  - `vscore`: The variance of the log-rank test statistic for the
    time-to-event endpoint.

  - `logRankStatistic`: The log-rank test Z-statistic for the
    time-to-event endpoint.

- `rawdataBIN` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for the binary endpoint
  for selected replications, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage under consideration.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `treatmentGroup`: The treatment group (1 or 2) for the subject.

  - `survivalTime`: The underlying survival time for the time-to-event
    endpoint for the subject.

  - `dropoutTime`: The underlying dropout time for the time-to-event
    endpoint for the subject.

  - `trtDiscTime`: The underlying treatment discontinuation time for the
    binary endpoint for the subject.

  - `trtDurUpperLimit`: The maximum protocol-specified treatment
    duration for the subject based on the treatment group assignment.

  - `ptfu1Time`:The underlying assessment time for the binary endpoint
    for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization for the binary endpoint for the subject.

  - `latentResponse`: The underlying latent response variable for the
    binary endpoint for the subject, which determines the response
    status for the binary endpoint at PTFU1 visit.

  - `responder`: Whether the subject is a responder for the binary
    endpoint.

  - `source`: The source of the determination of responder status for
    the binary endpoint: = 1 based on the underlying latent response
    variable, = 2 based on the occurrence of the time-to-event endpoint
    before the assessment time of the binary endpoint (imputed as a
    non-responder), = 3 based on the dropout before the assessment time
    of the binary endpoint (imputed as a non-responder), = 4 excluded
    from analysis due to administrative censoring.

- `rawdataTTE` (exists if `maxNumberOfRawDatasetsPerStage` is a positive
  integer): A data frame for subject-level data for the time-to-event
  endpoint for selected replications, containing the following
  variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage under consideration.

  - `analysisTime`: The time for the stage since trial start.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `treatmentGroup`: The treatment group (1 or 2) for the subject.

  - `survivalTime`: The underlying survival time for the time-to-event
    endpoint for the subject.

  - `dropoutTime`: The underlying dropout time for the time-to-event
    endpoint for the subject.

  - `timeUnderObservation`: The time under observation since
    randomization for the time-to-event endpoint for the subject.

  - `event`: Whether the subject experienced the event for the
    time-to-event endpoint.

  - `dropoutEvent`: Whether the subject dropped out for the
    time-to-event endpoint.

## Details

We consider dual primary endpoints with endpoint 1 being a binary
endpoint and endpoint 2 being a time-to-event endpoint. The analyses of
endpoint 1 will be based on calendar times, while the analyses of
endpoint 2 will be based on the number of events. Therefore, the
analyses of the two endpoints are not at the same time points. The
correlation between the two endpoints is characterized by the global
odds ratio of the Plackett copula. In addition, the time-to-event
endpoint will render the binary endpoint as a non-responder, and so does
the dropout. In addition, the treatment discontinuation will impact the
number of available subjects for analysis. The administrative censoring
will exclude subjects from the analysis of the binary endpoint.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
tcut <- c(0, 12, 36, 48)
surv <- c(1, 0.95, 0.82, 0.74)
lambda2 <- (log(surv[1:3]) - log(surv[2:4]))/(tcut[2:4] - tcut[1:3])

sim1 <- binary_tte_sim(
  kMax1 = 1,
  kMax2 = 2,
  accrualTime = seq(0, 8),
  accrualIntensity = 40/9 * seq(1, 9),
  piecewiseSurvivalTime = c(0, 12, 36),
  globalOddsRatio = 1,
  pi1 = 0.80,
  pi2 = 0.65,
  lambda1 = 0.65*lambda2,
  lambda2 = lambda2,
  gamma1 = -log(1-0.04)/12,
  gamma2 = -log(1-0.04)/12,
  delta1 = -log(1-0.02)/12,
  delta2 = -log(1-0.02)/12,
  upper1 = 15*28/30.4,
  upper2 = 12*28/30.4,
  n = 640,
  plannedTime = 20 + 15*28/30.4,
  plannedEvents = c(130, 173),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasetsPerStage = 1,
  seed = 314159,
  nthreads = 1)
```
