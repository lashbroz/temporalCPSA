# Build an ordered segment diagnostic matrix

Fit contiguous AD-TMP segmentations for a range of `k` values and return
a matrix whose rows are `k` values and whose columns are age-ordered
samples. This provides the data needed to draw a green age-ordered
diagnostic heatmap analogous to the manuscript age-class selection
display, while preserving contiguity for every fitted row.

## Usage

``` r
ageTMP_segment_diagnostic_matrix(
  tmp_matrix,
  sample_age,
  k_values = 2:6,
  min_n = 10,
  scale_rows = FALSE
)
```

## Arguments

- tmp_matrix:

  Numeric AD-TMP matrix with features in rows and samples in columns.

- sample_age:

  Numeric sample ages. If named, names are matched to
  `colnames(tmp_matrix)` after sample-ID normalization. If unnamed, the
  vector must be in the same order as `colnames(tmp_matrix)`.

- k_values:

  Integer vector of segment counts to fit.

- min_n:

  Minimum number of samples allowed in each segment.

- scale_rows:

  Whether to center and scale AD-TMP rows before scoring segments. Set
  this to `FALSE` when `tmp_matrix` has already been row-scaled.

## Value

A list with `matrix`, `ordered_samples`, and per-`k` `fits`.
