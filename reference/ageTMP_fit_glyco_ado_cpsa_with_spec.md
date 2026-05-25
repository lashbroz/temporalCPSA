# Fit glyco ADO CPSA models with an explicit specification

Fit glyco ADO CPSA models with an explicit specification

## Usage

``` r
ageTMP_fit_glyco_ado_cpsa_with_spec(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  spec = ageTMP_glyco_ado_cpsa_spec(cohort = "discovery"),
  n_cores = 1,
  progress = FALSE,
  progress_every = 500
)
```

## Arguments

- feature_matrix:

  Numeric glycopeptide-by-sample matrix.

- survival_data:

  Prepared survival data.

- features:

  Optional glycopeptide subset.

- spec:

  CPSA specification. Defaults to
  [`ageTMP_glyco_ado_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_glyco_ado_cpsa_spec.md).

- n_cores:

  Number of cores for fitting.

- progress:

  Whether to print fitting progress.

- progress_every:

  Progress chunk size.

## Value

A CPSA result data frame.
