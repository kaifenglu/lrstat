# Power and Sample Size for Cochran-Armitage Trend Test for Ordered Multi-Sample Binomial Response

Obtains the power given sample size or obtains the sample size given
power for the Cochran-Armitage trend test for ordered multi-sample
binomial response.

## Usage

``` r
getDesignOrderedBinom(
  beta = NA_real_,
  n = NA_real_,
  ngroups = NA_integer_,
  pi = NA_real_,
  w = NA_real_,
  allocationRatioPlanned = NA_integer_,
  rounding = TRUE,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- ngroups:

  The number of treatment groups.

- pi:

  The response probabilities for the treatment groups.

- w:

  The scores assigned to the treatment groups. This should reflect the
  ordinal nature of the treatment groups, e.g. dose levels. Defaults to
  equally spaced scores.

- allocationRatioPlanned:

  Allocation ratio for the treatment groups.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designOrderedBinom` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The maximum number of subjects.

- `ngroups`: The number of treatment groups.

- `pi`: The response probabilities for the treatment groups.

- `w`: The scores assigned to the treatment groups.

- `trendstat`: The Cochran-Armitage trend test statistic.

- `allocationRatioPlanned`: Allocation ratio for the treatment groups.

- `rounding`: Whether to round up sample size.

## Details

An ordered multi-sample binomial response design is used to test whether
the response probabilities differ across multiple treatment groups. The
null hypothesis is that the response probabilities are equal across all
treatment groups, while the alternative hypothesis is that the response
probabilities are ordered, i.e. the response probability increases with
the treatment group index. The Cochran-Armitage trend test is used to
test this hypothesis. This test effectively regresses the response
probabilities against treatment group scores, and test whether the slope
of the regression line is significantly different from zero.

The trend parameter is defined as \$\$\theta = \sum\_{g=1}^{G} r_g
(w_g - \bar{w}) \pi_g\$\$ where \\G\\ is the number of treatment groups,
\\r_g\\ is the randomization probability for treatment group \\g\\,
\\w_g\\ is the score assigned to treatment group \\g\\, \\\pi_g\\ is the
response probability for treatment group \\g\\, and \\\bar{w} =
\sum\_{g=1}^{G} r_g w_g\\ is the weighted average score across all
treatment groups.

Since \\\hat{\theta}\\ is a linear combination of the estimated response
probabilities, its variance is given by \$\$Var(\hat{\theta}) =
\frac{1}{n}\sum\_{g=1}^{G} r_g (w_g - \bar{w})^2 \pi_g(1-\pi_g)\$\$
where \\n\\ is the total sample size.

The sample size is chosen such that the power to reject the null
hypothesis is at least \\1-\beta\\ for a given significance level
\\\alpha\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignOrderedBinom(
  beta = 0.1, ngroups = 3, pi = c(0.1, 0.25, 0.5), alpha = 0.05))
#>   alpha     power  n ngroups trendstat
#> 1  0.05 0.9011121 75       3        10
```
