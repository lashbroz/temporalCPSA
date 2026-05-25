# Build manuscript STable4 discovery-cohort CPSA columns

This modality-neutral helper fits the cDisc/discovery CPSA model to a
public feature matrix and returns the signed log10 p-value columns used
in STable4.

## Usage

``` r
ageTMP_build_sa_discovery_cpsa(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho"),
  features = NULL,
  mode = c("manuscript", "standardized")
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- modality:

  Molecular modality accepted by
  [`ageTMP_load_feature_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_feature_matrix.md).

- features:

  Optional feature subset.

- mode:

  Reproducibility mode. `"manuscript"` preserves the sex-specific
  mutation missingness behavior used in the manuscript CPSA scripts.

## Value

A data frame with signed log10 p-value columns for male and female.
