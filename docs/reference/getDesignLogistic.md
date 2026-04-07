# Power and Sample Size for Logistic Regression

Obtains the power given sample size or obtains the sample size given
power for logistic regression of a binary response given the covariate
of interest and other covariates.

## Usage

``` r
getDesignLogistic(
  beta = NA_real_,
  n = NA_real_,
  ncovariates = NA_integer_,
  nconfigs = NA_integer_,
  x = NA_real_,
  pconfigs = NA_real_,
  corr = 0,
  oddsratios = NA_real_,
  responseprob = NA_real_,
  rounding = TRUE,
  alpha = 0.05
)
```

## Arguments

- beta:

  The type II error.

- n:

  The total sample size.

- ncovariates:

  The number of covariates.

- nconfigs:

  The number of configurations of discretized covariate values.

- x:

  The matrix of covariate values.

- pconfigs:

  The vector of probabilities for the configurations.

- corr:

  The multiple correlation between the predictor and other covariates.
  Defaults to 0.

- oddsratios:

  The odds ratios for one unit increase in the covariates.

- responseprob:

  The response probability in the full model when all predictor
  variables are equal to their means.

- rounding:

  Whether to round up sample size. Defaults to 1 for sample size
  rounding.

- alpha:

  The two-sided significance level. Defaults to 0.05.

## Value

An S3 class `designLogistic` object with the following components:

- `power`: The power to reject the null hypothesis.

- `alpha`: The two-sided significance level.

- `n`: The total sample size.

- `ncovariates`: The number of covariates.

- `nconfigs`: The number of configurations of discretized covariate
  values.

- `x`: The matrix of covariate values.

- `pconfigs`: The vector of probabilities for the configurations.

- `corr`: The multiple correlation between the predictor and other
  covariates.

- `oddsratios`: The odds ratios for one unit increase in the covariates.

- `responseprob`: The response probability in the full model when all
  predictor variables are equal to their means.

- `effectsize`: The effect size for the chi-square test.

- `rounding`: Whether to round up sample size.

## Details

We consider the logistic regression of a binary response variable \\Y\\
on a set of predictor variables \\x = (x_1,\ldots,x_K)^T\\ with \\x_1\\
being the covariate of interest: \\\log \frac{P(Y_i=1)}{1 - P(Y_i = 1)}
= \psi_0 + x_i^T \psi,\\ where \\\psi = (\psi_1,\ldots,\psi_K)^T\\.
Similar to Self et al (1992), we assume that all covariates are either
inherently discrete or discretized from continuous distributions (e.g.
using the quantiles). Let \\m\\ denote the total number of
configurations of the covariate values. Let \$\$\pi_i = P(x = x_i), i =
1,\ldots, m\$\$ denote the probabilities for the configurations of the
covariates under independence. The likelihood ratio test statistic for
testing \\H_0: \psi_1 = 0\\ can be approximated by a noncentral
chi-square distribution with one degree of freedom and noncentrality
parameter \$\$\Delta = 2 \sum\_{i=1}^m \pi_i \[b'(\theta_i)(\theta_i -
\theta_i^\*) - \\b(\theta_i) - b(\theta_i^\*)\\\],\$\$ where
\$\$\theta_i = \psi_0 + \sum\_{j=1}^{k} \psi_j x\_{ij},\$\$
\$\$\theta_i^\* = \psi_0^\* + \sum\_{j=2}^{k} \psi_j^\* x\_{ij},\$\$ for
\\\psi_0^\* = \psi_0 + \psi_1 \mu_1\\, and \\\psi_j^\* = \psi_j\\ for
\\j=2,\ldots,K\\. Here \\\mu_1\\ is the mean of \\x_1\\, e.g., \\\mu_1 =
\sum_i \pi_i x\_{i1}.\\ In addition, by formulating the logistic
regression in the framework of generalized linear models, \$\$b(\theta)
= \log(1 + \exp(\theta)),\$\$ and \$\$b'(\theta) =
\frac{\exp(\theta)}{1 + \exp(\theta)}.\$\$

The regression coefficients \\\psi\\ can be obtained by taking the log
of the odds ratios for the covariates. The intercept \\\psi_0\\ can be
derived as \$\$\psi_0 = \log(\bar{\mu}/(1- \bar{\mu})) - \sum\_{j=1}^{K}
\psi_j \mu_j,\$\$ where \\\bar{\mu}\\ denotes the response probability
when all predictor variables are equal to their means.

Finally, let \\\rho\\ denote the multiple correlation between the
predictor and other covariates. The noncentrality parameter of the
chi-square test is adjusted downward by multiplying by \\1-\rho^2\\.

## References

Steven G. Self, Robert H. Mauritsen and Jill Ohara. Power calculations
for likelihood ratio tests in generalized linear models. Biometrics
1992; 48:31-39.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# two ordinal covariates
x1 = c(5, 10, 15, 20)
px1 = c(0.2, 0.3, 0.3, 0.2)

x2 = c(2, 4, 6)
px2 = c(0.4, 0.4, 0.2)

# discretizing a normal distribution with mean 4 and standard deviation 2
nbins = 10
x3 = qnorm(((1:nbins) - 0.5)/nbins)*2 + 4
px3 = rep(1/nbins, nbins)

# combination of covariate values
nconfigs = length(x1)*length(x2)*length(x3)
x = expand.grid(x3 = x3, x2 = x2, x1 = x1)
x = as.matrix(x[, ncol(x):1])

# probabilities for the covariate configurations under independence
pconfigs = as.numeric(px1 %x% px2 %x% px3)

# convert the odds ratio for the predictor variable in 5-unit change
# to the odds ratio in 1-unit change
(design1 <- getDesignLogistic(
  beta = 0.1, ncovariates = 3,
  nconfigs = nconfigs,
  x = x,
  pconfigs = pconfigs,
  oddsratios = c(1.2^(1/5), 1.4, 1.3),
  responseprob = 0.25,
  alpha = 0.1))
#>   alpha     power    n ncovariates corr responseprob oddsratio  effectsize
#> 1   0.1 0.9001583 1369           3    0         0.25  1.037137 0.006259349
```
