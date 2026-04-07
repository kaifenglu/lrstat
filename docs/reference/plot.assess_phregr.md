# Plot method for assess_phregr objects

Generate line plots of observed and simulated paths of standardized
score processes of a Cox model.

## Usage

``` r
# S3 method for class 'assess_phregr'
plot(x, nsim = 20, ...)
```

## Arguments

- x:

  An object of class `assess_phregr`.

- nsim:

  The number of simulation samples used in the plot.

- ...:

  Ensures that all arguments starting from "..." are named.

## Value

A list of ggplot2 objects for the line plots, one for each covariate.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>
