# Build manuscript STable4 protein reference-cohort CPSA columns

Public-facing alias for
[`ageTMP_build_sa_protein_cdisc()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_sa_protein_cdisc.md).
It returns the manuscript `SA-Protein-cDisc-Ref` reference-cohort signed
log10 p-value columns, preserving the published column naming
convention.

## Usage

``` r
ageTMP_build_sa_protein_reference(
  data_dir = "data",
  genes = NULL,
  mode = c("manuscript", "standardized")
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- genes:

  Optional gene subset for testing.

- mode:

  Reproducibility mode. `"manuscript"` preserves the original
  sex-specific mutation missingness behavior used by the manuscript
  protein CPSA scripts: male mutation NAs are set to zero and female
  mutation NAs are preserved. `"standardized"` zero-fills missing
  mutation calls for both sexes.

## Value

A data frame with manuscript reference-cohort signed log10 p-value
columns for male and female.
