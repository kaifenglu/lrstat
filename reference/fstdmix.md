# Adjusted p-Values for Standard Mixture Gatekeeping Procedures

Obtains the adjusted p-values for the standard gatekeeping procedures
for multiplicity problems involving serial and parallel logical
restrictions.

## Usage

``` r
fstdmix(
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

Alex Dmitrienko and Ajit C Tamhane. Mixtures of multiple testing
procedures for gatekeeping applications in clinical trials. Statistics
in Medicine. 2011; 30(13):1473–1488.

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
fstdmix(p, family, serial, parallel, gamma, test = "hommel",
        exhaust = FALSE)
#> [1] 0.024250 0.013600 0.033875 0.024250 0.046250 0.024250 0.081400 0.033875
```
