# Combine TMP matrices for multi-omic sample clustering

Select dynamic rows within each molecular data type, restrict to common
samples, and stack the matrices for sample-level clustering. Row names
are prefixed by modality so that repeated feature identifiers remain
unique.

## Usage

``` r
ageTMP_combine_tmp_matrices(
  tmp_matrices,
  proportion = 0.5,
  max_features = 5000,
  common_samples = NULL,
  scale_rows = TRUE
)
```

## Arguments

- tmp_matrices:

  Named list of numeric feature-by-sample TMP matrices.

- proportion:

  Proportion of rows to keep after ranking by standard deviation.

- max_features:

  Maximum number of rows to keep.

- common_samples:

  Optional sample IDs to require. If `NULL`, the intersection of matrix
  column names is used.

- scale_rows:

  Whether to row-center and row-scale the selected matrix.

## Value

A numeric stacked TMP matrix.
