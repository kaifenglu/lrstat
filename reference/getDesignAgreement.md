# Power and Sample Size for Cohen's kappa

Obtains the power given sample size or obtains the sample size given
power for Cohen's kappa.

## Usage

``` r
getDesignAgreement(
  beta = NA_real_,
  n = NA_real_,
  ncats = NA_integer_,
  kappaH0 = NA_real_,
  kappa = NA_real_,
  p1 = NA_real_,
  p2 = NA_real_,
  rounding = TRUE,
  alpha = 0.025
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- ncats:

  The number of categories.

- kappaH0:

  The kappa coefficient under the null hypothesis.

- kappa:

  The kappa coefficient under the alternative hypothesis.

- p1:

  The marginal probabilities for the first rater.

- p2:

  The marginal probabilities for the second rater. Defaults to be equal
  to the marginal probabilities for the first rater if not provided.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The one-sided significance level. Defaults to 0.025.

## Value

An S3 class `designAgreement` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The one-sided significance level.

- `n`: The total sample size.

- `ncats`: The number of categories.

- `kappaH0`: The kappa coefficient under the null hypothesis.

- `kappa`: The kappa coefficient under the alternative hypothesis.

- `p1`: The marginal probabilities for the first rater.

- `p2`: The marginal probabilities for the second rater.

- `piH0`: The cell probabilities that maximize the variance of estimated
  kappa under H0.

- `pi`: The cell probabilities that maximize the variance of estimated
  kappa under H1.

- `rounding`: Whether to round up sample size.

## Details

The kappa coefficient is defined as \$\$\kappa = \frac{\pi_o -
\pi_e}{1 - \pi_e},\$\$ where \\\pi_o = \sum_i \pi\_{ii}\\ is the
observed agreement, and \\\pi_e = \sum_i \pi\_{i.} \pi\_{.i}\\ is the
expected agreement by chance.

By Fleiss et al. (1969), the variance of \\\hat{\kappa}\\ is given by
\$\$Var(\hat{\kappa}) = \frac{v_1}{n},\$\$ where \$\$v_1 = \frac{Q_1 +
Q_2 - Q3 - Q4}{(1-\pi_e)^4},\$\$ \$\$Q_1 = \pi_o(1-\pi_e)^2,\$\$ \$\$Q_2
= (1-\pi_o)^2 \sum_i \sum_j \pi\_{ij}(\pi\_{i.} + \pi\_{.j})^2,\$\$
\$\$Q_3 = 2(1-\pi_o)(1-\pi_e) \sum_i \pi\_{ii}(\pi\_{i.} +
\pi\_{.i}),\$\$ \$\$Q_4 = (\pi_o \pi_e - 2\pi_e + \pi_o)^2.\$\$

Given \\\kappa\\ and marginals \\\\(\pi\_{i.}, \pi\_{.i}):
i=1,\ldots,k\\\\, we obtain \\\pi_o\\. The only unknowns are the double
summation in \\Q_2\\ and the single summation in \\Q_3\\.

We find the optimal configuration of cell probabilities that yield the
maximum variance of \\\hat{\kappa}\\ by treating the problem as a linear
programming problem with constraints to match the given marginal
probabilities and the observed agreement and ensure that the cell
probabilities are nonnegative. This is an extension of Flack et al.
(1988) by allowing unequal marginal probabilities of the two raters.

We perform the optimization under both the null and alternative
hypotheses to obtain \\\max Var(\hat{\kappa} \| \kappa = \kappa_0)\\ and
\\\max Var(\hat{\kappa} \| \kappa = \kappa_1)\\ for a single subject,
and then calculate the sample size or power according to the following
equation: \$\$\sqrt{n} \|\kappa - \kappa_0\| = z\_{1-\alpha} \sqrt{\max
Var(\hat{\kappa} \| \kappa = \kappa_0)} + z\_{1-\beta} \sqrt{\max
Var(\hat{\kappa} \| \kappa = \kappa_1)}.\$\$

## References

V. F. Flack, A. A. Afifi, and P. A. Lachenbruch. Sample size
determinations for the two rater kappa statistic. Psychometrika 1988;
53:321-325.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignAgreement(
  beta = 0.2, n = NA, ncats = 4, kappaH0 = 0.4, kappa = 0.6,
  p1 = c(0.1, 0.2, 0.3, 0.4), p2 = c(0.15, 0.2, 0.24, 0.41),
  rounding = TRUE, alpha = 0.05))
#>   alpha   power  n kappaH0 kappa
#> 1  0.05 0.80064 82     0.4   0.6
```
