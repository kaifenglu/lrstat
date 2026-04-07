# Update Graph for Graphical Approaches

Updates the weights and transition matrix for graphical approaches.

## Usage

``` r
updateGraph(w, G, I, j)
```

## Arguments

- w:

  The current vector of weights for elementary hypotheses.

- G:

  The current transition matrix.

- I:

  The set of indices for yet to be rejected hypotheses.

- j:

  The hypothesis to remove from index set `I`.

## Value

A list containing the new vector of weights, the new transition matrix
for the graph, and the new set of indices of yet to be rejected
hypotheses.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
updateGraph(w = c(0.5, 0.5, 0, 0),
            G = matrix(c(0, 0.5, 0.5, 0,  0.5, 0, 0, 0.5,
                         0, 1, 0, 0,  1, 0, 0, 0),
                       nrow=4, ncol=4, byrow=TRUE),
            I = c(1, 2, 3, 4),
            j = 1)
#> $w
#> [1] 0.00 0.75 0.25 0.00
#> 
#> $G
#>      [,1] [,2]      [,3]      [,4]
#> [1,]    0  0.0 0.0000000 0.0000000
#> [2,]    0  0.0 0.3333333 0.6666667
#> [3,]    0  1.0 0.0000000 0.0000000
#> [4,]    0  0.5 0.5000000 0.0000000
#> 
#> $I
#> [1] 2 3 4
#> 
```
