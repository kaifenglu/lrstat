# Power and Sample Size for the Wilcoxon Test for Two-Sample Ordinal Response

Obtains the power given sample size or obtains the sample size given
power for the Wilcoxon test for two-sample ordinal response.

## Usage

``` r
getDesignTwoOrdinal(
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

  The number of categories of the ordinal response.

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

  The significance level. Defaults to 0.025.

## Value

An S3 class `designTwoOrdinal` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The maximum number of subjects.

- `ncats`: The number of categories of the ordinal response.

- `pi1`: The prevalence of each category for the treatment group.

- `pi2`: The prevalence of each category for the control group.

- `meanscore1`: The mean midrank score for the treatment group.

- `meanscore2`: The mean midrank score for the control group.

- `allocationRatioPlanned`: Allocation ratio for the active treatment
  versus control.

- `rounding`: Whether to round up sample size.

## Details

A two-sample ordinal response design is used to test whether the ordinal
response distributions differ between two treatment arms. Let
\\\pi\_{gi}\\ denote the prevalence of category \\i\\ in group \\g\\,
where \\g=1\\ represents the treatment group and \\g=2\\ represents the
control group.

The parameter of interest is \$\$\theta = \sum\_{i=1}^{C} w_i
(\pi\_{1i} - \pi\_{2i})\$\$ where \\w_i\\ is the midrank score for
category \\i\\. The Z-test statistic is given by \$\$Z =
\hat{\theta}/\sqrt{Var(\hat{\theta})}\$\$ where \\\hat{\theta}\\ is the
estimate of \\\theta\\.

The midrank score \\w_i\\ for category \\i\\ is calculated as: \$\$w_i =
\sum\_{j=1}^{i} \pi_j - 0.5\pi_i\$\$ where \\\pi_i = r\pi\_{1i} +
(1-r)\pi\_{2i}\\ denotes the average prevalence of category \\i\\ across
both groups, and \\r\\ is the randomization probability for the active
treatment.

To understand the midrank score, consider \\n\pi_i\\ subjects in
category \\i\\. The midrank score is the average rank of these subjects:
\$\$s_i = \frac{1}{n\pi_i} \sum\_{j=1}^{n\pi_i} ( n\pi_1 + \cdots +
n\pi\_{i-1} + j)\$\$ This simplifies to \$\$s_i = n\left(\sum\_{j=1}^{i}
\pi_j - 0.5\pi_i\right) + \frac{1}{2}\$\$ By dividing by \\n\\ and
ignoring \\\frac{1}{2n}\\, we obtain the midrank score \\w_i\\.

The variance of \\\hat{\theta}\\ can be derived from the multinomial
distributions and is given by \$\$Var(\hat{\theta}) =
\frac{1}{n}\sum\_{g=1}^{2} \frac{1}{r_g} \left\\\sum\_{i=1}^{C}
w_i^2\pi\_{gi} - \left(\sum\_{i=1}^{C} w_i\pi\_{gi}
\right)^2\right\\\$\$ where \\r_g\\ is the randomization probability for
group \\g\\.

The sample size is chosen such that the power to reject the null
hypothesis is at least \\1-\beta\\ for a given significance level
\\\alpha\\.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignTwoOrdinal(
  beta = 0.1, ncats = 4, pi1 = c(0.55, 0.3, 0.1),
  pi2 = c(0.214, 0.344, 0.251), alpha = 0.025))
#>   alpha     power  n ncats meanscore1 meanscore2
#> 1 0.025 0.9030195 67     4   26.40554   40.59446
```
