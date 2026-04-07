# Hedges' g Effect Size

Obtains Hedges' g estimate and confidence interval of effect size.

## Usage

``` r
hedgesg(tstat, m, ntilde, cilevel = 0.95)
```

## Arguments

- tstat:

  The value of the t-test statistic for comparing two treatment
  conditions.

- m:

  The degrees of freedom for the t-test.

- ntilde:

  The normalizing sample size to convert the standardized treatment
  difference to the t-test statistic, i.e.,
  `tstat = sqrt(ntilde)*meanDiff/stDev`.

- cilevel:

  The confidence interval level. Defaults to 0.95.

## Value

A data frame with the following variables:

- `tstat`: The value of the `t` test statistic.

- `m`: The degrees of freedom for the t-test.

- `ntilde`: The normalizing sample size to convert the standardized
  treatment difference to the t-test statistic.

- `g`: Hedges' `g` effect size estimate.

- `varg`: Variance of `g`.

- `lower`: The lower confidence limit for effect size.

- `upper`: The upper confidence limit for effect size.

- `cilevel`: The confidence interval level.

## Details

Hedges' \\g\\ is an effect size measure commonly used in meta-analysis
to quantify the difference between two groups. It's an improvement over
Cohen's \\d\\, particularly when dealing with small sample sizes.

The formula for Hedges' \\g\\ is \$\$g = c(m) d\$\$ where \\d\\ is
Cohen's \\d\\ effect size estimate, and \\c(m)\\ is the bias correction
factor, \$\$d = (\hat{\mu}\_1 - \hat{\mu}\_2)/\hat{\sigma}\$\$ \$\$c(m)
= 1 - \frac{3}{4m-1}.\$\$ Since \\c(m) \< 1\\, Cohen's \\d\\
overestimates the true effect size, \\\delta = (\mu_1 - \mu_2)/\sigma.\\
Since \$\$t = \sqrt{\tilde{n}} d\$\$ we have \$\$g =
\frac{c(m)}{\sqrt{\tilde{n}}} t\$\$ where \\t\\ has a noncentral \\t\\
distribution with \\m\\ degrees of freedom and noncentrality parameter
\\\sqrt{\tilde{n}} \delta\\.

The asymptotic variance of \\g\\ can be approximated by \$\$Var(g) =
\frac{1}{\tilde{n}} + \frac{g^2}{2m}.\$\$ The confidence interval for
\\\delta\\ can be constructed using normal approximation.

For two-sample mean difference with sample size \\n_1\\ for the
treatment group and \\n_2\\ for the control group, we have \\\tilde{n} =
\frac{n_1n_2}{n_1+n_2}\\ and \\m=n_1+n_2-2\\ for pooled variance
estimate.

## References

Larry V. Hedges. Distribution theory for Glass's estimator of effect
size and related estimators. Journal of Educational Statistics 1981;
6:107-128.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
n1 = 7
n2 = 8
meanDiff = 0.444
stDev = 1.201
m = n1+n2-2
ntilde = n1*n2/(n1+n2)
tstat = sqrt(ntilde)*meanDiff/stDev

hedgesg(tstat, m, ntilde)
#>       tstat  m   ntilde         g      varg      lower    upper cilevel
#> 1 0.7143127 13 3.733333 0.3479453 0.2725135 -0.6752113 1.371102    0.95
```
