# Fit manuscript glyco ADO-span CPSA models

The glyco SA scripts test the combined ADO feature effect using only PED
and ADO samples, with PED as the reference age class. This wrapper
preserves that model while reusing the package's modality-neutral CPSA
engine.

## Usage

``` r
ageTMP_fit_glyco_ado_cpsa(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  cohort = c("discovery", "reference"),
  scale_features = TRUE,
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

- cohort:

  One of `"discovery"` or `"reference"`.

- scale_features:

  Whether to row-scale feature values inside the model. The manuscript
  scripts scaled rows immediately before model fitting.

- n_cores:

  Number of cores for fitting.

- progress:

  Whether to print fitting progress.

- progress_every:

  Progress chunk size.

## Value

A CPSA result data frame.
