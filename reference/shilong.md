# The randomized clinical trial SHIVA data in long format from the ipcwswitch package

The original SHIdat data set contains an anonymized excerpt of data from
the SHIVA01 trial. This was the first randomized clinical trial that
aimed at comparing molecularly targeted therapy based on tumor profiling
(MTA) versus conventional therapy (CT) for advanced cancer. Patients
were randomly assigned to receive the active or control treatment and
may switch to the other arm or subsequent anti-cancer therapy upon
disease progression. The restructured data is in the long format.

- `id`:

  The patient's identifier

- `tstart`:

  The start of the time interval

- `tstop`:

  The end of the time interval

- `event`:

  Whether the patient died at the end of the interval

- `agerand`:

  The patient's age (in years) at randomization

- `sex.f`:

  The patients' gender, either Male or Female

- `tt_Lnum`:

  The number of previous lines of treatment

- `rmh_alea.c`:

  The Royal Marsden Hospital score segregated into two categories

- `pathway.f`:

  The molecular pathway altered (the hormone receptors pathway, the
  PI3K/ AKT/mTOR pathway, and the RAF/MEK pathway)

- `bras.f`:

  The patient's randomized arm, either MTA or CT

- `ps`:

  The ECOG performance status

- `ttc`:

  The presence of concomitant treatments

- `tran`:

  The use of platelet transfusions

- `dpd`:

  The relative day of a potential progression

- `dco`:

  The relative day of treatment switching

- `ady`:

  The relative day of the latest news

- `dcut`:

  The relative day of administrative cutoff

- `pd`:

  Whether the patient had disease progression

- `co`:

  Whether the patient switched treatment

## Usage

``` r
shilong
```

## Format

An object of class `data.frame` with 602 rows and 19 columns.
