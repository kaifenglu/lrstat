# lrstat

`lrstat` provides power and sample size methods for non-proportional
hazards and many other clinical trial designs.

The package is built around weighted log-rank methodology for
time-to-event group sequential designs, with flexible accrual,
event/dropout modeling, error-spending boundaries, and simulation
support. It also includes design and inference tools for continuous,
binary, count, and equivalence settings, including adaptive and
multi-arm/multi-stage extensions.

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("kaifenglu/lrstat")
```

## Key Features

- Weighted log-rank design and power under non-proportional hazards
  using the Fleming-Harrington class of weights.
- Group sequential design support with Lan-DeMets error spending and
  boundary calculations.
- Analytical and simulation-based operating characteristics for power,
  event counts, accrual duration, and follow-up duration.
- Support for adaptive group sequential updates (sample size, timing,
  error spending, future looks).
- Design modules for continuous, binary, count, and time-to-event
  endpoints, including equivalence settings.
- Specialized methods for crossover, repeated measures (MMRM), exact
  methods, and multi-arm/multi-stage workflows.
- Built-in datasets and a Shiny interface.

## Typical Workflow

1.  Specify design assumptions: accrual profile, event/dropout hazards,
    allocation, number/timing of looks, and spending function.
2.  Compute design characteristics with weighted log-rank functions such
    as
    [`lrpower()`](https://kaifenglu.github.io/lrstat/reference/lrpower.md),
    [`lrsamplesize()`](https://kaifenglu.github.io/lrstat/reference/lrsamplesize.md),
    [`getBound()`](https://kaifenglu.github.io/lrstat/reference/getBound.md),
    and related utilities.
3.  Optionally evaluate robustness via simulation (for example,
    [`lrsim()`](https://kaifenglu.github.io/lrstat/reference/lrsim.md))
    and compare alternative design scenarios.
4.  Summarize final design choices and explore sensitivity to delayed
    effects, accrual changes, or follow-up constraints.

## Minimal Time-to-Event Example

The example below computes power for a two-look group sequential trial
with a delayed treatment effect and FH(0,1) weighting.

``` r
library(lrstat)

fit <- lrpower(
  kMax = 2,
  informationRates = c(0.8, 1),
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  allocationRatioPlanned = 1,
  accrualTime = seq(0, 9),
  accrualIntensity = c(26 / 9 * seq(1, 9), 26),
  piecewiseSurvivalTime = c(0, 6),
  lambda1 = c(0.0533, 0.0309),
  lambda2 = c(0.0533, 0.0533),
  gamma1 = -log(1 - 0.05) / 12,
  gamma2 = -log(1 - 0.05) / 12,
  accrualDuration = 22,
  followupTime = 18,
  fixedFollowup = FALSE,
  rho1 = 0,
  rho2 = 1
)

fit
```

## Additional Design Families

`lrstat` includes broad design support beyond weighted log-rank
settings, including:

- Continuous endpoints (for example, mean difference/ratio and
  MMRM-based designs).
- Binary endpoints (including exact methods and interval-based
  dose-finding tools such as mTPI-2 and BOIN helpers).
- Count endpoints (including negative binomial settings with censoring).
- Equivalence and adaptive group sequential variants across supported
  endpoint types.

See the reference index for the full function catalog.

## Shiny App

Launch the interactive application:

``` r
library(lrstat)
runShinyApp_lrstat()
```

## Documentation

- Website: <https://kaifenglu.github.io/lrstat/>
- Function reference: <https://kaifenglu.github.io/lrstat/reference/>
- Vignettes: <https://kaifenglu.github.io/lrstat/articles/>
- Source code: <https://github.com/kaifenglu/lrstat>
- Issue tracker: <https://github.com/kaifenglu/lrstat/issues>

## Citation

If you use `lrstat` in analyses, reports, or publications, please cite
the package and relevant methodological references documented in the
function help pages and vignettes.
