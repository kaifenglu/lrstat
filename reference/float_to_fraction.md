# Converting a decimal to a fraction

Converts a decimal to a fraction based on the algorithm from
http://stackoverflow.com/a/5128558/221955.

## Usage

``` r
float_to_fraction(x, tol = 1e-06)
```

## Arguments

- x:

  The fraction in decimal form.

- tol:

  The tolerance level for the conversion error.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
float_to_fraction(5/3)
#> [1] 5 3
```
