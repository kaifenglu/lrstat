# Power and Sample Size for Unordered Multi-Sample Binomial Response

Obtains the power given sample size or obtains the sample size given
power for the chi-square test for unordered multi-sample binomial
response.

## Usage

``` r
getDesignUnorderedBinom(
  beta = NA_real_,
  n = NA_real_,
  ngroups = NA_integer_,
  pi = NA_real_,
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

- allocationRatioPlanned:

  Allocation ratio for the treatment groups.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designUnorderedBinom` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The maximum number of subjects.

- `ngroups`: The number of treatment groups.

- `pi`: The response probabilities for the treatment groups.

- `effectsize`: The effect size for the chi-square test.

- `allocationRatioPlanned`: Allocation ratio for the treatment groups.

- `rounding`: Whether to round up sample size.

## Details

A multi-sample binomial response design is used to test whether the
response probabilities differ among multiple treatment arms. Let
\\\pi\_{g}\\ denote the response probability in group \\g =
1,\ldots,G\\, where \\G\\ is the total number of treatment groups.

The chi-square test statistic is given by \$\$X^2 = \sum\_{g=1}^{G}
\sum\_{i=1}^{2} \frac{(n\_{gi} - n\_{g+}n\_{+i}/n)^2}{n\_{g+}
n\_{+i}/n}\$\$ where \\n\_{gi}\\ is the number of subjects in category
\\i\\ for group \\g\\, \\n\_{g+}\\ is the total number of subjects in
group \\g\\, and \\n\_{+i}\\ is the total number of subjects in category
\\i\\ across all groups, and \\n\\ is the total sample size.

Let \\r_g\\ denote the randomization probability for group \\g\\, and
define the weighted average response probability across all groups as
\$\$\bar{\pi} = \sum\_{g=1}^{G} r_g \pi_g\$\$

- Under the null hypothesis, \\X^2\\ follows a chi-square distribution
  with \\G-1\\ degrees of freedom.

- Under the alternative hypothesis, \\X^2\\ follows a non-central
  chi-square distribution with non-centrality parameter \$\$\lambda = n
  \sum\_{g=1}^{G} \frac{r_g (\pi\_{g} - \bar{\pi})^2} {\bar{\pi}
  (1-\bar{\pi})}\$\$

The sample size is chosen such that the power to reject the null
hypothesis is at least \\1-\beta\\ for a given significance level
\\\alpha\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignUnorderedBinom(
  beta = 0.1, ngroups = 3, pi = c(0.1, 0.25, 0.5), alpha = 0.05))
#>   alpha     power  n ngroups effectsize
#> 1  0.05 0.9019526 95       3  0.1340629
```
