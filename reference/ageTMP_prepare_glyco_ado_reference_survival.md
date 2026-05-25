# Prepare reference glyco ADO-span survival data

Prepare reference glyco ADO-span survival data

## Usage

``` r
ageTMP_prepare_glyco_ado_reference_survival(
  clinical,
  mutation = NULL,
  glyco_sample_ids = NULL,
  sex = c("Male", "Female"),
  max_day = 2000,
  require_complete_mutation = TRUE
)
```

## Arguments

- clinical:

  Clinical data frame, usually from
  [`ageTMP_load_reference_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_reference_clinical.md).

- mutation:

  Optional sample-by-gene mutation indicator matrix.

- glyco_sample_ids:

  Glyco trajectory sample IDs available for modeling.

- sex:

  One of `"Male"` or `"Female"`.

- max_day:

  Administrative survival truncation in days.

- require_complete_mutation:

  For public reference clinical tables that contain `has.complete.mut`,
  restrict to samples with complete mutation covariates. This matches
  the manuscript reference-cohort CPSA model frame.

## Value

A prepared reference survival data frame.
