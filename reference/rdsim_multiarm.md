# Risk Difference Simulation for Multi-Arm Multi-Stage Design

Simulate a multi-arm multi-stage design for binary responses using
risk-difference Wald statistics with closed testing.

## Usage

``` r
rdsim_multiarm(
  M = 2,
  kMax = 1,
  criticalValues = NA,
  futilityBounds = NULL,
  riskDiffH0s = 0,
  allocations = 1,
  pis = NULL,
  nullVariance = TRUE,
  n = NA,
  plannedSubjects = NA,
  maxNumberOfIterations = 1000,
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

  Numeric matrix of dimension \\kMax \times M\\ giving the by-look
  critical values for the closed testing procedure. The first column is
  used for the level-M test and the last column for the level-1 test.

- futilityBounds:

  Numeric vector of length \\kMax - 1\\ giving the futility boundaries
  on the Wald-statistic scale for the first \\kMax - 1\\ looks. At an
  interim look, the study stops for futility if all active treatment
  arms fall below the futility boundary. If omitted, no interim futility
  stopping is applied.

- riskDiffH0s:

  Scalar or numeric vector of length \\M\\. Risk differences under
  \\H_0\\ for each active arm versus the common control. Defaults to 0.

- allocations:

  Integer or integer vector of length \\M + 1\\. Number of subjects per
  arm within a randomization block. A single value implies equal
  allocation; defaults to 1. The first \\M\\ elements refer to the
  active arms and the last element refers to the common control.

- pis:

  Numeric vector of length \\M + 1\\. Each element corresponds to the
  response rate for a treatment arm. The first \\M\\ elements refer to
  the active arms and the last element refers to the common control.

- nullVariance:

  Whether to use the variance under the null or the empirical variance
  under the alternative.

- n:

  Planned total sample size across all active arms and control.

- plannedSubjects:

  Numeric vector of length \\kMax\\ giving the planned cumulative sample
  size at each look for the first active arm and the common control
  combined.

- maxNumberOfIterations:

  Number of Monte Carlo replications. Defaults to 1000.

- seed:

  Random seed for reproducibility.

- nthreads:

  Number of threads for parallel simulation. Use 0 to accept the default
  RcppParallel behavior.

## Value

An S3 object of class `"rdsim_multiarm"` with these components:

- `overview`: A list summarizing trial-level results and settings:

  - `overallReject`: Overall probability of rejecting the null by trial
    end.

  - `overallFutility`: Overall probability of stopping for futility by
    trial end.

  - `rejectPerStage`: Probability of rejecting the null for each active
    arm at each stage.

  - `futilityPerStage`: Probability of futility stopping for each active
    arm at each stage.

  - `cumulativeRejection`: Cumulative probability of rejection by stage.

  - `cumulativeFutility`: Cumulative futility stopping probability by
    stage.

  - `numberOfEvents`: Cumulative event counts by stage and arm.

  - `numberOfSubjects`: Cumulative enrollments by stage and arm.

  - `expectedNumberOfEvents`: Expected cumulative events at trial end.

  - `expectedNumberOfSubjects`: Expected cumulative enrollments at trial
    end.

  - `criticalValues`: The input matrix of by-level critical values.

  - `futilityBounds`: The input futility boundaries for each stage.

  - `riskDiffH0s`: The input risk differences under \\H_0\\.

  - `nullVariance`: Whether to use variance under \\H_0\\.

  - `numberOfIterations`: Number of simulation iterations performed.

  - `n`: Planned total sample size.

  - `allocations`: The input allocation ratios.

  - `responseRates`: The input response rates for each arm.

  - `plannedSubjects`: The input planned cumulative sample size at each
    look for the first active arm and the common control combined.

  - `M`: Number of active arms.

  - `kMax`: Number of sequential looks.

- `sumdata1`: Data frame summarizing each iteration, stage, and
  treatment group:

  - `iterationNumber`, `stopStage`, `stageNumber`, `treatmentGroup`,
    `accruals`, `events`, `phat`.

  - For each stage the final row summarizes the overall study (all arms
    combined).

- `sumdata2`: Data frame summarizing test statistics by iteration,
  stage, and active arm:

  - `iterationNumber`, `stopStage`, `stageNumber`, `activeArm`,
    `totalAccruals`, `totalEvents`, `riskDiff`, `vriskDiff`,
    `riskDiffZ`, `reject`, `futility`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(sim1 <- rdsim_multiarm(
  M = 2,
  kMax = 3,
  criticalValues = matrix(c(3.880, 2.747, 2.275,
                            3.710, 2.511, 1.993), 3, 2),
  futilityBounds =  c(0.043, 1.194),
  pis = c(0.25, 0.30, 0.20),
  n = 486,
  plannedSubjects = c(146, 292, 324),
  maxNumberOfIterations = 10000,
  seed = 314159,
  nthreads = 0))
#>                                                  
#> Multi-arm multi-stage design for risk difference 
#> Overall power: 0.4586                            
#> Number of active arms: 2                         
#> Number of looks: 3                               
#> n: 486, null variance: TRUE                      
#> Number of simulations: 10000                     
#>                                                  
#> By level critical boundaries
#>         Level 2 Level 1
#> Stage 1   3.880   3.710
#> Stage 2   2.747   2.511
#> Stage 3   2.275   1.993
#> 
#> Cumulative probability of rejection or futility by treatment
#>         Reject Active 1 Reject Active 2 Overall Rejection Futility
#> Stage 1          0.0011          0.0055            0.0064   0.0551
#> Stage 2          0.0556          0.2193            0.2328   0.1833
#> Stage 3          0.1186          0.4326            0.4586   0.5414
#> 
#> Detailed probability of trial termination at each look
#>         Reject 1 Active Reject 2 Actives Overall Rejection Futility Continue
#> Stage 1          0.0062           0.0002            0.0064   0.0551   0.9385
#> Stage 2          0.1845           0.0419            0.2264   0.1282   0.5839
#> Stage 3          0.1753           0.0505            0.2258   0.3581   0.0000
#> Total            0.3660           0.0926            0.4586   0.5414       NA
#> 
#> Overall probability of rejection by set of active arms
#>   Set of active arms Probability of rejection
#> 1               none                   0.5414
#> 2                  1                   0.0260
#> 3                  2                   0.3400
#> 4                1,2                   0.0926
#> 
#>                            Active 1 Active 2 Control Total
#> Expected # events              37.7     45.2    30.1 113.0
#> Expected # subjects           150.9    150.9   150.9 452.6
#> 
#>                            Active 1 Active 2 Control Total
#> Number of events   Stage 1     18.3     21.9    14.6  54.7
#> Number of events   Stage 2     36.6     43.9    28.9 109.4
#> Number of events   Stage 3     40.6     48.0    32.5 121.1
#> Number of subjects Stage 1     73.0     73.0    73.0 219.0
#> Number of subjects Stage 2    146.0    146.0   146.0 438.0
#> Number of subjects Stage 3    162.0    162.0   162.0 486.0
```
