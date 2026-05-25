# Collapse matrix rows by feature identifier

Collapse matrix rows by feature identifier

## Usage

``` r
ageTMP_collapse_matrix_by_feature(mat, feature, method = "mean")
```

## Arguments

- mat:

  Numeric feature-by-sample matrix.

- feature:

  Character vector assigning each row of `mat` to a collapsed feature
  identifier.

- method:

  Summary method. Currently only `"mean"` is supported.

## Value

A numeric matrix with one row per unique feature.
