# The binary data from Cox and Snell (1989, pp. 10-11).

The dataset consits of the number of ingots not ready for rolling and
the number of ingots ready for rolling for a number of combinations of
heating time and soaking time.

## Usage

``` r
ingots
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 25
rows and 4 columns.

## Details

- `Heat`:

  The heating time

- `Soak`:

  The soaking time

- `NotReady`:

  Response indicator, with a value 1 for units not ready for rolling
  (event) and a value of 0 for units ready for rolling (nonevent)

- `Freq`:

  The frequency of occurrence of each combination of `Heat`, `Soak`, and
  `NotReady`
