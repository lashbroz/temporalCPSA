# Prepare discovery glyco ADO-span survival data

This wraps the existing discovery survival preparation with manuscript
glyco defaults: mutation missingness is zero-filled, age classes are
restricted to PED/ADO, and samples are filtered to a supplied
glycopeptide sample set.

## Usage

``` r
ageTMP_prepare_glyco_ado_discovery_survival(
  clinical,
  mutation = NULL,
  glyco_sample_ids = NULL,
  sex = c("Male", "Female"),
  max_day = 2000
)
```

## Arguments

- clinical:

  Discovery clinical data.

- mutation:

  Optional sample-by-gene mutation indicator matrix.

- glyco_sample_ids:

  Glyco sample IDs available for modeling.

- sex:

  One of `"Male"` or `"Female"`.

- max_day:

  Administrative survival truncation in days.

## Value

A prepared survival data frame.
