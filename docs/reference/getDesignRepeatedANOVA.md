# Power and Sample Size for Repeated-Measures ANOVA

Obtains the power and sample size for one-way repeated measures analysis
of variance. Each subject takes all treatments in the longitudinal
study.

## Usage

``` r
getDesignRepeatedANOVA(
  beta = NA_real_,
  n = NA_real_,
  ngroups = 2,
  means = NA_real_,
  stDev = 1,
  corr = 0,
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

- means:

  The treatment group means.

- stDev:

  The total standard deviation.

- corr:

  The correlation among the repeated measures.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designRepeatedANOVA` object with the following components:

- `power`: The power to reject the null hypothesis that there is no
  difference among the treatment groups.

- `alpha`: The two-sided significance level.

- `n`: The number of subjects.

- `ngroups`: The number of treatment groups.

- `means`: The treatment group means.

- `stDev`: The total standard deviation.

- `corr`: The correlation among the repeated measures.

- `effectsize`: The effect size.

- `rounding`: Whether to round up sample size.

## Details

Let \\y\_{ij}\\ denote the measurement under treatment condition \\j
(j=1,\ldots,k)\\ for subject \\i (i=1,\ldots,n)\\. Then \$\$y\_{ij} =
\alpha + \beta_j + b_i + e\_{ij}\$\$ where \\b_i\\ denotes the subject
random effect, \\b_i \sim N(0, \sigma_b^2)\\ and \\e\_{ij} \sim N(0,
\sigma_e^2)\\ denotes the within-subject residual. If we set \\\beta_k =
0\\, then \\\alpha\\ is the mean of the last treatment (control), and
\\\beta_j\\ is the difference in means between the \\j\\th treatment and
the control for \\j=1,\ldots,k-1\\.

The repeated measures have a compound symmetry covariance structure. Let
\\\sigma^2 = \sigma_b^2 + \sigma_e^2\\, and \\\rho =
\frac{\sigma_b^2}{\sigma_b^2 + \sigma_e^2}\\. Then \\Var(y_i) = \sigma^2
\\(1-\rho) I_k + \rho 1_k 1_k^T\\\\. Let \\X_i\\ denote the design
matrix for subject \\i\\. Let \\\theta = (\alpha, \beta_1, \ldots,
\beta\_{k-1})^T\\. It follows that \$\$Var(\hat{\theta}) =
\left(\sum\_{i=1}^{n} X_i^T V_i^{-1} X_i\right)^{-1}.\$\$ It can be
shown that \$\$Var(\hat{\beta}) = \frac{\sigma^2 (1-\rho)}{n}
(I\_{k-1} + 1\_{k-1} 1\_{k-1}^T).\$\$ It follows that \\\hat{\beta}^T
\hat{V}\_{\hat{\beta}}^{-1} \hat{\beta} \sim F\_{k-1,(n-1)(k-1),
\lambda}\\ where the noncentrality parameter for the \\F\\ distribution
is \$\$\lambda = \beta^T V\_{\hat{\beta}}^{-1} \beta = \frac{n
\sum\_{j=1}^{k} (\mu_j - \bar{\mu})^2}{\sigma^2(1-\rho)}.\$\$

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
(design1 <- getDesignRepeatedANOVA(
  beta = 0.1, ngroups = 4, means = c(1.5, 2.5, 2, 0),
  stDev = 5, corr = 0.2, alpha = 0.05))
#>   alpha     power  n ngroups stDev corr effectsize
#> 1  0.05 0.9027338 83       4     5  0.2      0.175
```
