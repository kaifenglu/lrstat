# The liver data used in SAS PROC PHREG documentation examples.

This data set contains information on 418 patients with primary biliary
cirrhosis.

## Usage

``` r
liver
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 418
rows and 7 columns.

## Details

- `Time`:

  The follow-up time in years from the time of registration. The event
  could be liver transplantation, death, or the end of the study,
  whichever came first

- `Status`:

  A censoring indicator, where a value of 1 indicating a death event,
  and 0 indicating a censored observation (the patient survived past the
  observation time)

- `Age`:

  The patient's age in years

- `Albumin`:

  Serum albumin level in g/dl

- `Bilirubin`:

  Serum bilirubin level in mg/dl

- `Edema`:

  Edema status, where a value of 0 indicates no edema, 0.5 indicates
  edema successfully treated with diuretics, and 1 indicates edema
  despite diuretic therapy

- `Protime`:

  Prothrombin time in seconds
