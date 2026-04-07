# Simulated CONCORDE trial data from the rpsftm package

Patients were randomly assigned to receive treatment immediately or
deferred, and those in the deferred arm could cross over and receive
treatment. The primary endpoint was time to disease progression.

- `id`:

  Patient identification number

- `def`:

  Indicator that the participant was assigned to the deferred treatment
  arm

- `imm`:

  Indicator that the participant was assigned to the immediate treatment
  arm

- `censyrs`:

  The censoring time, in years, corresponding to the close of study
  minus the time of entry for each patient

- `xo`:

  Indicator that crossover occurred

- `xoyrs`:

  The time, in years, from entry to switching, or 0 for patients in the
  immediate arm

- `prog`:

  Indicator of disease progression (1), or censoring (0)

- `progyrs`:

  Time, in years, from entry to disease progression or censoring

- `entry`:

  The time of entry into the study, measured in years from the date of
  randomisation

## Usage

``` r
immdef
```

## Format

An object of class `data.frame` with 1000 rows and 9 columns.
