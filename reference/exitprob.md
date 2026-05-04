# Stagewise Exit Probabilities

Obtains the stagewise exit probabilities for both efficacy and futility
stopping.

## Usage

``` r
exitprob(b, a = NA_real_, theta = 0L, I = NA_real_)
```

## Arguments

- b:

  Upper boundaries on the z-test statistic scale.

- a:

  Lower boundaries on the z-test statistic scale. Defaults to
  `c(rep(-8.0, kMax-1), b[kMax])` if left unspecified, where
  `kMax = length(b)`.

- theta:

  Stagewise parameter of interest, e.g., `-U/V` for weighted log-rank
  test, where `U` is the mean and `V` is the variance of the weighted
  log-rank test score statistic at each stage. For proportional hazards
  and conventional log-rank test, use the scalar input,
  `theta = -log(HR)`. Defaults to 0 corresponding to the null
  hypothesis.

- I:

  Stagewise cumulative information, e.g., `V`, the variance of the
  weighted log-rank test score statistic at each stage. For conventional
  log-rank test, information can be approximated by `phi*(1-phi)*D`,
  where `phi` is the probability of being allocated to the active arm,
  and `D` is the total number of events at each stage. Defaults to
  `seq(1, kMax)` if left unspecified.

## Value

A list of stagewise exit probabilities:

- `exitProbUpper`: The vector of efficacy stopping probabilities

- `exitProbLower`: The vector of futility stopping probabilities.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
exitprob(b = c(3.471, 2.454, 2.004), theta = -log(0.6),
         I = c(50, 100, 150)/4)
#> $exitProbUpper
#> [1] 0.0479605 0.4927446 0.3327286
#> 
#> $exitProbLower
#> [1] 5.302375e-23 2.428197e-26 1.265665e-01
#> 

exitprob(b = c(2.963, 2.359, 2.014),
         a = c(-0.264, 0.599, 2.014),
         theta = c(0.141, 0.204, 0.289),
         I = c(81, 121, 160))
#> $exitProbUpper
#> [1] 0.04513264 0.40913028 0.44590828
#> 
#> $exitProbLower
#> [1] 0.06263793 0.02066973 0.01652114
#> 
```
