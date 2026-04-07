# Analysis of Simon's Bayesian Basket Trials

Obtains the prior and posterior probabilities for Simon's Bayesian
basket discovery trials.

## Usage

``` r
simonBayesAnalysis(
  nstrata = NA_integer_,
  r = NA_real_,
  n = NA_real_,
  lambda = NA_real_,
  gamma = NA_real_,
  phi = NA_real_,
  plo = NA_real_
)
```

## Arguments

- nstrata:

  The number of strata.

- r:

  The vector of number of responders across strata.

- n:

  The vector of number of subjects across strata.

- lambda:

  The prior probability that the drug activity is homogeneous across
  strata.

- gamma:

  The prior probability that the drug is active in a stratum.

- phi:

  The response probability for an active drug.

- plo:

  The response probability for an inactive drug.

## Value

A list containing the following five components:

- `case`: The matrix with each row corresponding to a combination of
  drug activity over strata represented by the columns.

- `prior_case`: The vector of joint prior probabilities for the
  stratum-specific response rates.

- `prior_stratum`: The vector of marginal prior probabilities for the
  stratum-specific response rates.

- `post_case`: The vector of joint posterior probabilities for the
  stratum-specific response rates.

- `post_stratum`: The vector of marginal posterior probabilities for the
  stratum-specific response rates.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
a = simonBayesAnalysis(
  nstrata = 10,
  r = c(8,0,1,1,6,2,0,0,3,3),
  n = c(19,10,26,8,14,7,8,5,4,14),
  lambda = 0.5, gamma = 0.33,
  phi = 0.35, plo = 0.15)

a$post_stratum
#>  [1] 0.950120630 0.032389925 0.001458153 0.148349241 0.895742213 0.408975723
#>  [7] 0.054091888 0.113259675 0.820611914 0.244634490
```
