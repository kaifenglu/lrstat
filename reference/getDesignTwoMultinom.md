# Power and Sample Size for Difference in Two-Sample Multinomial Responses

Obtains the power given sample size or obtains the sample size given
power for difference in two-sample multinomial responses.

## Usage

``` r
getDesignTwoMultinom(
  beta = NA_real_,
  n = NA_real_,
  ncats = NA_integer_,
  pi1 = NA_real_,
  pi2 = NA_real_,
  allocationRatioPlanned = 1,
  rounding = TRUE,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- ncats:

  The number of categories of the multinomial response.

- pi1:

  The prevalence of each category for the treatment group. Only need to
  specify the valued for the first `ncats-1` categories.

- pi2:

  The prevalence of each category for the control group. Only need to
  specify the valued for the first `ncats-1` categories.

- allocationRatioPlanned:

  Allocation ratio for the active treatment versus control. Defaults to
  1 for equal randomization.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designTwoMultinom` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The maximum number of subjects.

- `ncats`: The number of categories of the multinomial response.

- `pi1`: The prevalence of each category for the treatment group.

- `pi2`: The prevalence of each category for the control group.

- `effectsize`: The effect size for the chi-square test.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `rounding`: Whether to round up sample size.

## Details

A two-arm multinomial response design is used to test whether the
prevalence of each category differs between two treatment arms. Let
\\\pi\_{gi}\\ denote the prevalence of category \\i\\ in group \\g\\,
where \\g=1\\ for the treatment group and \\g=2\\ for the control group.
The chi-square test statistic is given by \$\$X^2 = \sum\_{g=1}^{2}
\sum\_{i=1}^{C} \frac{(n\_{gi} - n\_{g+} n\_{+i}/n)^2}{n\_{g+}
n\_{+i}/n}\$\$ where \\n\_{gi}\\ is the number of subjects in category
\\i\\ for group \\g\\, \\n\_{g+}\\ is the total number of subjects in
group \\g\\, and \\n\_{+i}\\ is the total number of subjects in category
\\i\\ across both groups, and \\n\\ is the total sample size.

- Under the null hypothesis, \\X^2\\ follows a chi-square distribution
  with \\C-1\\ degrees of freedom.

- Under the alternative hypothesis, \\X^2\\ follows a non-central
  chi-square distribution with non-centrality parameter \$\$\lambda = n
  r (1-r) \sum\_{i=1}^{C} \frac{(\pi\_{1i} - \pi\_{2i})^2} {r
  \pi\_{1i} + (1-r)\pi\_{2i}}\$\$ where \\r\\ is the randomization
  probability for the active treatment.

The sample size is chosen such that the power to reject the null
hypothesis is at least \\1-\beta\\ for a given significance level
\\\alpha\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignTwoMultinom(
  beta = 0.1, ncats = 3, pi1 = c(0.3, 0.35),
  pi2 = c(0.2, 0.3), alpha = 0.05))
#>   alpha     power   n ncats effectsize
#> 1  0.05 0.9000174 503     3 0.02515837
```
