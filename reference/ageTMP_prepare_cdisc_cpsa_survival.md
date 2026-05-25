# Prepare reference-cohort clinical covariates for CPSA survival modeling

Prepare reference-cohort clinical covariates for CPSA survival modeling

## Usage

``` r
ageTMP_prepare_cdisc_cpsa_survival(
  clinical,
  mutation = NULL,
  protein_sample_ids = NULL,
  sex = NULL,
  age_classes = c("PED", "ADO", "YA", "ADULT"),
  max_day = 2000,
  mutation_na = c("zero", "preserve"),
  require_complete_mutation = TRUE
)
```

## Arguments

- clinical:

  Clinical data frame, usually from
  [`ageTMP_load_reference_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_reference_clinical.md).

- mutation:

  Optional sample-by-gene mutation indicator matrix.

- protein_sample_ids:

  Optional vector of sample IDs available in the protein matrix.

- sex:

  Optional sex filter, `"Male"` or `"Female"`.

- age_classes:

  Age classes to keep.

- max_day:

  Administrative survival truncation in days.

- mutation_na:

  How to handle missing mutation covariates after matching mutation
  calls to clinical samples. Use `"zero"` to treat missing calls as wild
  type/absent, or `"preserve"` to keep missing values as `NA`. The
  manuscript protein CPSA workflow used sex-specific behavior: missing
  mutation covariates were zero-filled in the male run and preserved in
  the female run.

- require_complete_mutation:

  For public reference clinical tables that contain `has.complete.mut`,
  restrict to samples with complete mutation covariates. This matches
  the manuscript reference-cohort CPSA model frame.

## Value

A data frame with survival outcome and model covariates.
