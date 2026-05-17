# Group Sequential Trials Using Bonferroni-Based Graphical Approaches

Obtains the test results for group sequential trials using graphical
approaches based on weighted Bonferroni tests.

## Usage

``` r
fseqbon(
  w,
  G,
  alpha = 0.025,
  kMax,
  typeAlphaSpending = NULL,
  parameterAlphaSpending = NULL,
  maxInformation = NULL,
  incidenceMatrix = NULL,
  k1,
  p,
  information,
  spendingTime = NULL,
  nthreads = 0
)
```

## Arguments

- w:

  The vector of initial weights for elementary hypotheses.

- G:

  The initial transition matrix.

- alpha:

  The significance level. Defaults to 0.025.

- kMax:

  The maximum number of stages.

- typeAlphaSpending:

  The vector of alpha spending functions for the hypotheses. Each
  element is one of the following: `"sfOF"` for O'Brien-Fleming type
  spending function, `"sfP"` for Pocock type spending function, `"sfKD"`
  for Kim & DeMets spending function, `"sfHSD"` for Hwang, Shi & DeCani
  spending function. Defaults to `"sfOF"` if not provided.

- parameterAlphaSpending:

  The vector of parameter values for the alpha spending functions for
  the hypotheses. Each element corresponds to the value of \\\rho\\ for
  `"sfKD"` or \\\gamma\\ for `"sfHSD"`. Defaults to missing if not
  provided.

- maxInformation:

  The vector of target maximum information for each hypothesis. Defaults
  to a vector of 1s if not provided.

- incidenceMatrix:

  The `kMax x m` incidence matrix indicating whether the specific
  hypothesis will be tested at the given look. Here \\m\\ is the number
  of hypotheses. If not provided, defaults to testing each hypothesis at
  all study looks.

- k1:

  The number of study looks at the interim analysis.

- p:

  The \\B k_1 \times m\\ matrix of raw p-values for each hypothesis by
  study look. Here \\B\\ is the number of replications for simulations.
  In other words, the number of rows must be a multiple of `k1`.

- information:

  The \\B k_1 \times m\\ matrix of observed information for each
  hypothesis by study look.

- spendingTime:

  The \\B k_1 \times m\\ matrix of spending time for alpha spending by
  study look. Each element must be between 0 and 1, and the spending
  time must be increasing by study look for each hypothesis. If not
  provided, it is the same as `informationRates` calculated from
  `information` and `maxInformation`.

- nthreads:

  The number of threads to use in simulations (0 means the default
  RcppParallel behavior).

## Value

A vector to indicate the first look the specific hypothesis is rejected
(0 if the hypothesis is not rejected).

## Details

The procedure allows the user to retest an unrejected hypothesis at
earlier looks if its weight increased after rejection of other
hypotheses. In this case, the procedure will return the first look at
which the specific hypothesis is rejected. If the hypothesis is not
rejected at any look, it will return 0.

## References

Willi Maurer and Frank Bretz. Multiple testing in group sequential
trials using graphical approaches. Statistics in Biopharmaceutical
Research. 2013; 5:311-320.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Case study from Maurer & Bretz (2013)

fseqbon(
  w = c(0.5, 0.5, 0, 0),
  G = matrix(c(0, 0.5, 0.5, 0,  0.5, 0, 0, 0.5,
               0, 1, 0, 0,  1, 0, 0, 0),
             nrow=4, ncol=4, byrow=TRUE),
  alpha = 0.025,
  kMax = 3,
  typeAlphaSpending = rep("sfOF", 4),
  maxInformation = rep(1, 4),
  k1 = 2,
  p = matrix(c(0.0062, 0.017, 0.009, 0.13,
               0.0002, 0.0035, 0.002, 0.06),
             nrow=2, ncol=4, byrow=TRUE),
  information = matrix(c(rep(1/3, 4), rep(2/3, 4)),
                       nrow=2, ncol=4, byrow=TRUE),
  nthreads = 1)
#> [1] 2 2 2 0

```
