# Adjusted p-Values for Modified Mixture Gatekeeping Procedures

Obtains the adjusted p-values for the modified gatekeeping procedures
for multiplicity problems involving serial and parallel logical
restrictions.

## Usage

``` r
fmodmix(
  p,
  family = NULL,
  serial,
  parallel = NULL,
  gamma,
  test = "hommel",
  exhaust = TRUE
)
```

## Arguments

- p:

  The raw p-values for elementary hypotheses.

- family:

  The matrix of family indicators for the hypotheses.

- serial:

  The matrix of serial rejection set for the hypotheses.

- parallel:

  The matrix of parallel rejection set for the hypotheses.

- gamma:

  The truncation parameters for each family. The truncation parameter
  for the last family is automatically set to 1.

- test:

  The component multiple testing procedure. Options include "holm",
  "hochberg", or "hommel". Defaults to "hommel".

- exhaust:

  Whether to use alpha-exhausting component testing procedure for the
  last family with active hypotheses. It defaults to `TRUE`.

## Value

A matrix of adjusted p-values.

## References

Alex Dmitrienko, George Kordzakhia, and Thomas Brechenmacher.
Mixture-based gatekeeping procedures for multiplicity problems with
multiple sequences of hypotheses. Journal of Biopharmaceutical
Statistics. 2016; 26(4):758–780.

George Kordzakhia, Thomas Brechenmacher, Eiji Ishida, Alex Dmitrienko,
Winston Wenxiang Zheng, and David Fuyuan Li. An enhanced mixture method
for constructing gatekeeping procedures in clinical trials. Journal of
Biopharmaceutical Statistics. 2018; 28(1):113–128.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
p <- c(0.0194, 0.0068, 0.0271, 0.0088, 0.0370, 0.0018, 0.0814, 0.0066)
family <- matrix(c(1, 1, 0, 0, 0, 0, 0, 0,
                   0, 0, 1, 1, 0, 0, 0, 0,
                   0, 0, 0, 0, 1, 1, 0, 0,
                   0, 0, 0, 0, 0, 0, 1, 1),
                 nrow=4, byrow=TRUE)

serial <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0,
                   1, 0, 0, 0, 0, 0, 0, 0,
                   0, 1, 0, 0, 0, 0, 0, 0,
                   0, 0, 1, 0, 0, 0, 0, 0,
                   0, 0, 0, 1, 0, 0, 0, 0,
                   0, 0, 0, 0, 1, 0, 0, 0,
                   0, 0, 0, 0, 0, 1, 0, 0),
                 nrow=8, byrow=TRUE)

parallel <- matrix(0, 8, 8)
gamma <- c(0.6, 0.6, 0.6, 1)
fmodmix(p, family, serial, parallel, gamma, test = "hommel", exhaust = TRUE)
#> [1] 0.02425 0.01360 0.03300 0.02425 0.03700 0.02425 0.08140 0.03300
```
