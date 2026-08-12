# Risk Difference Simulation for Phase 2/3 Seamless Design

Simulate phase 2/3 seamless design testing for risk difference.

## Usage

``` r
rdsim_seamless(
  M = 2,
  K = 1,
  rankp0 = 1,
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

  Number of active treatment arms in Phase 2.

- K:

  Number of sequential looks in Phase 3.

- rankp0:

  Integer rank in phase 2 used to select the active arm. `rankp0 = 1`
  selects the most efficacious arm by risk-difference statistic,
  `rankp0 = 2` selects the second most efficacious arm, and so on.

- criticalValues:

  Numeric vector of length \\K + 1\\ giving the critical value for the
  Wald statistic at each look (Look 1 through Look \\K + 1\\). Decision
  rule:

  - At Look 1, compute the Wald statistic for each active arm versus the
    common control. If the `rankp0`-th largest test statistic exceeds
    the Look 1 critical value, stop for efficacy.

  - If the Look 1 stopping rule is not met, select the active arm with
    the `rankp0`-th largest Wald statistic and continue with that arm
    only versus control at subsequent looks.

  - For each look \\j = 2,\ldots,K+1\\, compare the selected arm to
    control; if its Wald statistic exceeds the Look \\j\\ critical
    value, stop for efficacy; otherwise continue.

  - If no critical value is exceeded by Look \\K + 1\\, the procedure
    ends without rejection.

- futilityBounds:

  Numeric vector of length \\K\\ giving the futility boundaries for
  Phase 2 and the first \\K-1\\ looks in Phase 3. The study stops for
  futility:

  - in Phase 2 if the selected treatment arm crosses the phase-2
    futility boundary;

  - in Phase 3 if the selected arm crosses the futility boundary at an
    interim look; If omitted, no interim futility stopping is applied.

- riskDiffH0s:

  Numeric vector of length \\M\\. Risk differences under \\H_0\\ for
  each active arm versus the common control. Defaults to 0 for
  superiority tests.

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

  Numeric vector of length \\K + 1\\ giving the planned cumulative
  sample size at each look. Each entry refers to the combined number of
  patients for the first active arm and the common control.

- maxNumberOfIterations:

  Number of Monte Carlo replications. Defaults to 1000.

- seed:

  Random seed for reproducibility.

- nthreads:

  Number of threads for parallel simulation. Use 0 to accept the default
  RcppParallel behavior.

## Value

An S3 object of class `"rdsim_seamless"` with these components:

- `overview`: A list summarizing trial-level results and settings:

  - `selectionProb`: Probability of selecting each active arm at the
    prespecified rank at the end of phase 2.

  - `selectToStage2`: Probability of selecting each active arm to enter
    stage 2.

  - `selectAnyToStage2`: Probability of selecting any active arm to
    enter stage 2.

  - `rejectPerStage`: Probability of rejecting the null for each active
    arm at each stage.

  - `futilityPerStage`: Probability of futility stopping for each active
    arm at each stage.

  - `cumulativeRejection`: Cumulative probability of rejection by stage.

  - `cumulativeFutility`: Cumulative futility stopping probabilities by
    stage.

  - `numberOfEvents`: Cumulative event counts by stage, including events
    from all arms in stage 1 and events from the selected arm and
    control in later stages.

  - `numberOfSubjects`: Cumulative enrollments by stage. replications
    that reached that stage.

  - `overallReject`: Overall probability of rejecting the null by trial
    end.

  - `overallFutility`: Overall probability of stopping for futility by
    trial end.

  - `expectedNumberOfEvents`: Expected cumulative events at trial end.

  - `expectedNumberOfSubjects`: Expected cumulative enrollments at trial
    end.

  - `criticalValues`: The input critical values for each stage.

  - `futilityBounds`: The input futility boundaries for each stage.

  - `riskDiffH0s`: The input risk differences under \\H_0\\.

  - `nullVariance`: Whether to use variance under \\H_0\\.

  - `numberOfIterations`: Number of simulation iterations performed.

  - `n`: Planned total sample size.

  - `allocations`: The input allocation ratios.

  - `responseRates`: The input response rates for each arm.

  - `plannedSubjects`: The input planned cumulative sample size at each
    look for the first active arm and the common control combined.

  - `M`: Number of active arms in Phase 2.

  - `K`: Number of sequential looks in Phase 3.

  - `rankp0`: Prespecified rank used for phase-2 arm selection.

- `sumdata1`: Data frame summarizing each iteration, stage, and
  treatment group:

  - `iterationNumber`, `stopStage`, `stageNumber`, `treatmentGroup`,
    `accruals`, `events`, `phat`.

  - For each stage the final row summarizes the overall study (all arms
    combined).

- `summdata2`: Data frame summarizing test statistics by iteration,
  stage, and active arm:

  - `iterationNumber`, `selectedArm`, `stopStage`, `stageNumber`,
    `activeArm`, `totalAccruals`, `totalEvents`, `riskDiff`,
    `vriskDiff`, `riskDiffZ`, `reject`, `futility`.

  - For each active arm, total accruals and events refer to the combined
    counts for that arm and the common control at that stage.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(sim1 <- rdsim_seamless(
  M = 3,
  K = 1,
  criticalValues = c(8, 2.349),
  futilityBounds = 1.036,
  pis = c(0.22, 0.25, 0.35, 0.20),
  n = 1120,
  plannedSubjects = c(280, 560),
  maxNumberOfIterations = 10000,
  seed = 314159,
  nthreads = 0))
#>                                                                
#> Phase 2/3 seamless group-sequential design for risk difference 
#> Empirical power: 0.9192                                        
#> Number of active arms in phase 2: 3                            
#> Number of looks in phase 3: 1                                  
#> Selected rank in phase 2: 1                                    
#> Expected # events: 216.9                                       
#> Expected # subjects: 832.3                                     
#> n: 1120                                                        
#> Variance under H0: TRUE                                        
#> Number of simulations: 10000                                   
#>                                                                
#>                           Stage 1 Stage 2
#> Efficacy bounds (z-scale) 8.000   2.349  
#> Futility bounds (z-scale) 1.036   2.349  
#> 
#>                                           Arm 1  Arm 2  Arm 3
#> Selected at prespecified rank in phase 2 0.0081 0.0412 0.9507
#> 
#>                            Arm 1  Arm 2  Arm 3
#> Selected to enter stage 2 0.0074 0.0368 0.9284
#> 
#>                                          Probability
#> Any active arm selected to enter stage 2      0.9726
#> 
#>         Reject Active 1 Reject Active 2 Reject Active 3 Overall Rejection
#> Stage 1          0.0000          0.0000          0.0000            0.0000
#> Stage 2          0.0022          0.0177          0.8993            0.9192
#> Total            0.0022          0.0177          0.8993            0.9192
#>         Futility Continue
#> Stage 1   0.0274   0.9726
#> Stage 2   0.0534   0.0000
#> Total     0.0808       NA
#> 
#>  activeArm stage cumReject cumFutility nEvents nSubjects
#>          1     1    0.0000      0.0274   142.8     560.0
#>          1     2    0.0022      0.0326   200.9     840.0
#>          2     1    0.0000      0.0274   142.9     560.0
#>          2     2    0.0177      0.0465   205.9     840.0
#>          3     1    0.0000      0.0274   142.7     560.0
#>          3     2    0.8993      0.0565   219.7     840.0
#>    Overall     1    0.0000      0.0274   142.7     560.0
#>    Overall     2    0.9192      0.0808   219.0     840.0
```
