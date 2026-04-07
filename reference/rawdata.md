# A simulated time-to-event data set with 10 replications

A simulated data set with stratification and delayed treatment effect:

- `iterationNumber`:

  The iteration number

- `arrivalTime`:

  The enrollment time for the subject

- `stratum`:

  The stratum for the subject

- `treatmentGroup`:

  The treatment group for the subject

- `timeUnderObservation`:

  The time under observation since randomization

- `event`:

  Whether the subject experienced the event

- `dropoutEvent`:

  Whether the subject dropped out

## Usage

``` r
rawdata
```

## Format

An object of class `data.frame` with 4910 rows and 7 columns.
