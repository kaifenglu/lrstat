# Simulation of Simon's Bayesian Basket Trials

Obtains the simulated raw and summary data for Simon's Bayesian basket
discovery trials.

## Usage

``` r
simonBayesSim(
  p = NA_real_,
  accrualTime = 0L,
  accrualIntensity = NA_real_,
  stratumFraction = 1L,
  lambda = NA_real_,
  gamma = NA_real_,
  phi = NA_real_,
  plo = NA_real_,
  T = NA_real_,
  maxSubjects = NA_integer_,
  plannedSubjects = NA_integer_,
  maxNumberOfIterations = 1000L,
  maxNumberOfRawDatasets = 1L,
  seed = 0L
)
```

## Arguments

- p:

  The vector of true response probabilities across strata.

- accrualTime:

  A vector that specifies the starting time of piecewise Poisson
  enrollment time intervals. Must start with 0, e.g., `c(0, 3)` breaks
  the time axis into 2 accrual intervals: \\\[0, 3)\\ and \\\[3,
  \infty)\\.

- accrualIntensity:

  A vector of accrual intensities. One for each accrual time interval.

- stratumFraction:

  A vector of stratum fractions that sum to 1. Defaults to 1 for no
  stratification.

- lambda:

  The prior probability that the drug activity is homogeneous across
  strata.

- gamma:

  The prior probability that the drug is active in a stratum.

- phi:

  The response probability for an active drug.

- plo:

  The response probability for an inactive drug.

- T:

  The threshold for a conclusive posterior probability to stop
  enrollment.

- maxSubjects:

  The maximum total sample size.

- plannedSubjects:

  The planned cumulative number of subjects at each stage.

- maxNumberOfIterations:

  The number of simulation iterations. Defaults to 1000.

- maxNumberOfRawDatasets:

  The number of raw datasets to extract.

- seed:

  The seed to reproduce the simulation results. The seed from the
  environment will be used if left unspecified,

## Value

A list containing the following four components:

- `rawdata`: A data frame for subject-level data, containing the
  following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage number.

  - `subjectId`: The subject ID.

  - `arrivalTime`: The enrollment time for the subject.

  - `stratum`: The stratum for the subject.

  - `y`: Whether the subject was a responder (1) or nonresponder (0).

- `sumdata1`: A data frame for simulation and stratum-level summary
  data, containing the following variables:

  - `iterationNumber`: The iteration number.

  - `stageNumber`: The stage number.

  - `stratum`: The stratum number.

  - `active`: Whether the drug is active in the stratum.

  - `n`: The number of subjects in the stratum.

  - `r`: The number of responders in the stratum.

  - `posterior`: The posterior probability that the drug is active in
    the stratum.

  - `open`: Whether the stratum is still open for enrollment.

  - `positive`: Whether the stratum has been determined to be a positive
    stratum.

  - `negative`: Whether the stratum has been determined to be a negative
    stratum.

- `sumdata2`: A data frame for the simulation level summary data,
  containing the following variables:

  - `iterationNumber`: The iteration number.

  - `numberOfStrata`: The total number of strata.

  - `n_active_strata`: The number of active strata.

  - `true_positive`: The number of true positive strata.

  - `false_negative`: The number of false negative strata.

  - `false_positive`: The number of false positive strata.

  - `true_negative`: The number of true negative strata.

  - `n_indet_strata`: The number of indeterminate strata.

  - `numberOfSubjects`: The number of subjects.

- `overview`: A data frame for the summary across simulations,
  containing the following variables:

  - `numberOfStrata`: The total number of strata.

  - `n_active_strata`: The average number of active strata.

  - `true_positive`: The average number of true positive strata.

  - `false_negative`: The average number of false negative strata.

  - `false_positive`: The average number of false positive strata.

  - `true_negative`: The average number of true negative strata.

  - `n_indet_strata`: The average number of indeterminate strata.

  - `numberOfSubjects`: The average number of subjects.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
sim1 = simonBayesSim(
  p = c(0.25, 0.25, 0.05),
  accrualIntensity = 5,
  stratumFraction = c(1/3, 1/3, 1/3),
  lambda = 0.33, gamma = 0.5,
  phi = 0.25, plo = 0.05,
  T = 0.8, maxSubjects = 50,
  plannedSubjects = seq(5, 50, 5),
  maxNumberOfIterations = 1000,
  maxNumberOfRawDatasets = 1,
  seed = 314159)

sim1$overview
#>   numberOfStrata n_active_strata true_positive false_negative false_positive
#> 1              3               2         1.642          0.344          0.181
#>   true_negative n_indet_strata numberOfSubjects
#> 1         0.805          0.028           27.015
```
