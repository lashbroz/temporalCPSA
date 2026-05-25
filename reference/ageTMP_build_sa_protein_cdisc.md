# Build manuscript reference-cohort protein CPSA columns for STable4 comparison

Build manuscript reference-cohort protein CPSA columns for STable4
comparison

## Usage

``` r
ageTMP_build_sa_protein_cdisc(
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

A data frame with manuscript reference-cohort signed log10p columns for
male and female. Published `STable4` column names retain the `cdisc`
label because that is the manuscript source-table convention.
