# Power and Sample Size for One-Sample Multinomial Response

Obtains the power given sample size or obtains the sample size given
power for one-sample multinomial response.

## Usage

``` r
getDesignOneMultinom(
  beta = NA_real_,
  n = NA_real_,
  ncats = NA_integer_,
  piH0 = NA_real_,
  pi = NA_real_,
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

- piH0:

  The prevalence of each category under the null hypothesis. Only need
  to provide the values for the first `ncats-1` categories.

- pi:

  The prevalence of each category. Only need to provide the values for
  the first `ncats-1` categories.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designOneMultinom` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The maximum number of subjects.

- `ncats`: The number of categories of the multinomial response.

- `piH0`: The prevalence of each category under the null hypothesis.

- `pi`: The prevalence of each category.

- `effectsize`: The effect size for the chi-square test.

- `rounding`: Whether to round up sample size.

## Details

A single-arm multinomial response design is used to test whether the
prevalence of each category is different from the null hypothesis
prevalence. The null hypothesis is that the prevalence of each category
is equal to \\\pi\_{0i}\\, while the alternative hypothesis is that the
prevalence of each category is equal to \\\pi_i\\, for \\i=1,\ldots,C\\,
where \\C\\ is the number of categories.

The sample size is calculated based on the chi-square test for
multinomial response. The test statistic is given by \$\$X^2 =
\sum\_{i=1}^{C} \frac{(n_i - n\pi\_{0i})^2}{n\pi\_{0i}}\$\$ where
\\n_i\\ is the number of subjects in category \\i\\, and \\n\\ is the
total sample size.

- Under the null hypothesis, \\X^2\\ follows a chi-square distribution
  with \\C-1\\ degrees of freedom.

- Under the alternative hypothesis, \\X^2\\ follows a non-central
  chi-square distribution with non-centrality parameter \$\$\lambda = n
  \sum\_{i=1}^{C} \frac{(\pi_i - \pi\_{0i})^2}{\pi\_{0i}}\$\$

The sample size is chosen such that the power to reject the null
hypothesis is at least \\1-\beta\\ for a given significance level
\\\alpha\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignOneMultinom(
  beta = 0.1, ncats = 3, piH0 = c(0.25, 0.25),
  pi = c(0.3, 0.4), alpha = 0.05))
#>   alpha     power  n ncats effectsize
#> 1  0.05 0.9029864 71     3       0.18
```
