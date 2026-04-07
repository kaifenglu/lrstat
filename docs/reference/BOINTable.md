# BOIN Decision Table for Dose-Finding Trials

Generates the decision table for the Bayesian Optimal Interval (BOIN)
design, a widely used approach for dose-escalation trials that guides
dose-finding decisions based on observed toxicity rates.

## Usage

``` r
BOINTable(
  nMax = NA_integer_,
  pT = 0.3,
  phi1 = 0.6 * pT,
  phi2 = 1.4 * pT,
  a = 1,
  b = 1,
  pExcessTox = 0.95
)
```

## Arguments

- nMax:

  The maximum number of subjects allowed in a dose cohort.

- pT:

  The target toxicity probability. Defaults to 0.3.

- phi1:

  The lower equivalence limit for the target toxicity probability.

- phi2:

  The upper equivalence limit for the target toxicity probability.

- a:

  The prior toxicity shape parameter for the Beta prior.

- b:

  The prior non-toxicity shape parameter for the Beta prior.

- pExcessTox:

  The threshold for excessive toxicity. If the posterior probability
  that the true toxicity rate exceeds `pT` is greater than `pExcessTox`,
  the current and all higher doses will be excluded from further use to
  protect future participants. Defaults to 0.95.

## Value

An S3 class `BOINTable` object with the following components:

- `settings`: The input settings data frame with the following
  variables:

  - `nMax`: The maximum number of subjects in a dose cohort.

  - `pT`: The target toxicity probability.

  - `phi1`: The lower equivalence limit for target toxicity probability.

  - `phi2`: The upper equivalence limit for target toxicity probability.

  - `lambda1`: The lower decision boundary for observed toxicity
    probability.

  - `lambda2`: The upper decision boundary for observed toxicity
    probability.

  - `a`: The prior toxicity parameter for the beta prior.

  - `b`: The prior non-toxicity parameter for the beta prior.

  - `pExcessTox`: The threshold for excessive toxicity.

- `decisionDataFrame`: A data frame listing dose-finding decisions for
  each combination of sample size (`n`) and number of observed
  toxicities (`y`):

  - `n`: Cohort size.

  - `y`: Number of observed toxicities.

  - `decision`: Recommended action: escalate, de-escalate, or stay at
    the current dose.

- `decisionMatrix`: A matrix version of the decision table showing the
  recommended action based on the number of toxicities for each possible
  cohort size.

## References

Liu, S., & Yuan, Y. (2015). Bayesian optimal interval designs for phase
I clinical trials. Journal of the Royal Statistical Society: Series C
(Applied Statistics), 64(3), 507-523.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
BOINTable(nMax = 18, pT = 0.3, phi = 0.6*0.3, phi2 = 1.4*0.3)
#>                                        
#> Trial monitoring table for BOIN design 
#>                                        
#>    1 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
#> 0  E E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E
#> 1  D D  S  S  E  E  E  E  E  E  E  E  E  E  E  E  E  E
#> 2    D  D  D  D  S  S  S  E  E  E  E  E  E  E  E  E  E
#> 3      DU DU  D  D  D  D  S  S  S  S  E  E  E  E  E  E
#> 4         DU DU DU  D  D  D  D  D  S  S  S  S  S  E  E
#> 5            DU DU DU DU DU  D  D  D  D  S  S  S  S  S
#> 6               DU DU DU DU DU DU  D  D  D  D  D  S  S
#> 7                  DU DU DU DU DU DU DU  D  D  D  D  D
#> 8                     DU DU DU DU DU DU DU DU DU  D  D
#> 9                        DU DU DU DU DU DU DU DU DU DU
#> 10                          DU DU DU DU DU DU DU DU DU
#> 11                             DU DU DU DU DU DU DU DU
#> 12                                DU DU DU DU DU DU DU
#> 13                                   DU DU DU DU DU DU
#> 14                                      DU DU DU DU DU
#> 15                                         DU DU DU DU
#> 16                                            DU DU DU
#> 17                                               DU DU
#> 18                                                  DU
#>                                                                              
#> Rows represent number of toxicities                                          
#> Columns represent number of patients treated at current dose                 
#> E = Escalate to the next higher dose                                         
#> S = Stay at the current dose                                                 
#> D = De-escalate to the next lower dose                                       
#> DU = The current dose is unacceptably toxic                                  
#> Target toxicity: 0.3, phi1: 0.18, phi2: 0.42, lambda1: 0.236, lambda2: 0.359 
```
