# mTPI-2 Decision Table

Obtains the decision table for the modified toxicity probability
interval-2 (mTPI-2) design.

## Usage

``` r
mTPI2Table(
  nMax = NA_integer_,
  pT = 0.3,
  epsilon1 = 0.05,
  epsilon2 = 0.05,
  a = 1,
  b = 1,
  pExcessTox = 0.95
)
```

## Arguments

- nMax:

  The maximum number of subjects in a dose cohort.

- pT:

  The target toxicity probability. Defaults to 0.3.

- epsilon1:

  The lower equivalence margin from the target. Defaults to 0.05.

- epsilon2:

  The upper equivalence margin from the target. Defaults to 0.05.

- a:

  The prior toxicity parameter for the beta prior.

- b:

  The prior non-toxicity parameter for the beta prior.

- pExcessTox:

  The threshold for excessive toxicity, i.e., if Prob(p \> pT \| Data)
  \> pExcessTox, then the current and all higher doses will be excluded
  and never be used again in the remainder of the trial to avoid any
  other subjects receiving treatment at those doses. Defaults to 0.95.

## Value

An S3 class `mTPI2Table` object with the following components:

- `settings`: The input settings data frame with the following
  variables:

  - `nMax`: The maximum number of subjects in a dose cohort.

  - `pT`: The target toxicity probability.

  - `epsilon1`: The lower equivalence margin from the target.

  - `epsilon2`: The upper equivalence margin from the target.

  - `a`: The prior toxicity parameter for the beta prior.

  - `b`: The prior non-toxicity parameter for the beta prior.

  - `pExcessTox`: The threshold for excessive toxicity.

- `subintervals`: The subintervals of equal length in the mTPI-2 design.
  It includes the following variables:

  - `lower`: The lower bound of the subinterval.

  - `upper`: The upper bound of the subinterval.

  - `decision`: The dosing decision for the subinterval.

- `decisionDataFrame`: The decision data frame for the mTPI-2 design. It
  includes the following variables:

  - `n`: The sample size.

  - `y`: The number of toxicities.

  - `decision`: The dosing decision.

- `decisionMatrix`: The decision matrix corresponding to the decision
  data frame.

## References

Guo, W., Wang, S. J., Yang, S., Lynn, H., & Ji, Y. (2017). A Bayesian
interval dose-finding design addressing Ockham's razor: mTPI-2.
Contemporary Clinical Trials, 58, 23-33.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
mTPI2Table(nMax = 18, pT = 0.3, epsilon1 = 0.05, epsilon2 = 0.05)
#>                                          
#> Trial monitoring table for mTPI-2 design 
#>                                          
#>    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
#> 0  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E  E
#> 1  D  D  S  S  E  E  E  E  E  E  E  E  E  E  E  E  E  E
#> 2    DU  D  D  D  S  S  S  E  E  E  E  E  E  E  E  E  E
#> 3       DU DU  D  D  D  D  S  S  S  S  E  E  E  E  E  E
#> 4          DU DU DU  D  D  D  D  D  S  S  S  S  S  E  E
#> 5             DU DU DU DU DU  D  D  D  D  D  S  S  S  S
#> 6                DU DU DU DU DU DU  D  D  D  D  D  D  S
#> 7                   DU DU DU DU DU DU DU  D  D  D  D  D
#> 8                      DU DU DU DU DU DU DU DU DU  D  D
#> 9                         DU DU DU DU DU DU DU DU DU DU
#> 10                           DU DU DU DU DU DU DU DU DU
#> 11                              DU DU DU DU DU DU DU DU
#> 12                                 DU DU DU DU DU DU DU
#> 13                                    DU DU DU DU DU DU
#> 14                                       DU DU DU DU DU
#> 15                                          DU DU DU DU
#> 16                                             DU DU DU
#> 17                                                DU DU
#> 18                                                   DU
#>                                                              
#> Rows represent number of toxicities                          
#> Columns represent number of patients treated at current dose 
#> E = Escalate to the next higher dose                         
#> S = Stay at the current dose                                 
#> D = De-escalate to the next lower dose                       
#> DU = The current dose is unacceptably toxic                  
#> Target toxicity: 0.3, epsilon1: 0.05, epsilon2: 0.05         
```
