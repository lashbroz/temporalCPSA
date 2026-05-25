# Select dynamic TMP features

Select rows from a temporal molecular profile matrix by row-wise
standard deviation, then optionally row-center and row-scale the
selected matrix.

## Usage

``` r
ageTMP_select_dynamic_tmp_features(
  tmp_matrix,
  proportion = 0.5,
  max_features = 5000,
  scale_rows = TRUE
)
```

## Arguments

- tmp_matrix:

  Numeric feature-by-sample TMP matrix.

- proportion:

  Proportion of rows to keep after ranking by standard deviation.

- max_features:

  Maximum number of rows to keep.

- scale_rows:

  Whether to row-center and row-scale the selected matrix.

## Value

A numeric matrix containing selected dynamic TMP rows.
