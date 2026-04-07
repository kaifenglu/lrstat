# Simon's Two-Stage Design

Obtains Simon's two-stage minimax, admissible, and optimal designs.

## Usage

``` r
simon2stage(
  alpha = NA_real_,
  beta = NA_real_,
  piH0 = NA_real_,
  pi = NA_real_,
  n_max = 110L
)
```

## Arguments

- alpha:

  Type I error rate (one-sided).

- beta:

  Type II error rate (1-power).

- piH0:

  Response probability under the null hypothesis.

- pi:

  Response probability under the alternative hypothesis.

- n_max:

  Upper limit for sample size, defaults to 110.

## Value

A data frame containing the following variables:

- `piH0`: Response probability under the null hypothesis.

- `pi`: Response probability under the alternative hypothesis.

- `alpha`: The specified one-sided significance level.

- `beta`: The specified type II error.

- `n`: Total sample size.

- `n1`: Stage 1 sample size.

- `r1`: Futility boundary for stage 1.

- `r`: Futility boundary for stage 2.

- `EN0`: Expected sample size under the null hypothesis.

- `attainedAlpha`: Attained type 1 error.

- `power`: Attained power.

- `PET0`: Probability of early stopping under the null hypothesis.

- `w_lower`: Lower bound of the interval for `w`.

- `w_upper`: Upper bound of the interval for `w`.

- `design`: Description of the design, e.g., minimax, admissible, or
  optimal.

Here `w` is the weight in the objective function: `w*n + (1-w)*EN0`.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
simon2stage(0.05, 0.2, 0.1, 0.3)
#>   piH0  pi alpha beta  n n1 r1 r      EN0 attainedAlpha attainedPower      PET0
#> 1  0.1 0.3  0.05  0.2 25 15  1 5 19.50957    0.03280867     0.8017006 0.5490430
#> 2  0.1 0.3  0.05  0.2 26 12  1 5 16.77397    0.03596715     0.8047804 0.6590023
#> 3  0.1 0.3  0.05  0.2 27 11  1 5 15.84229    0.03950052     0.8061954 0.6973569
#> 4  0.1 0.3  0.05  0.2 29 10  1 5 15.01412    0.04708631     0.8050629 0.7360989
#>     w_lower   w_upper     design
#> 1 0.7323055 1.0000000    Minimax
#> 2 0.4823155 0.7323055 Admissible
#> 3 0.2928288 0.4823155 Admissible
#> 4 0.0000000 0.2928288    Optimal
```
