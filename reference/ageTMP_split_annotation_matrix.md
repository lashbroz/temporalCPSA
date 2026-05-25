# Split a molecular table into annotation and numeric matrix components

Split a molecular table into annotation and numeric matrix components

## Usage

``` r
ageTMP_split_annotation_matrix(
  data,
  annotation_cols,
  row_id = NULL,
  normalize_colnames = TRUE
)
```

## Arguments

- data:

  Molecular data frame with feature annotation columns followed by
  sample columns.

- annotation_cols:

  Integer or character vector identifying annotation columns.

- row_id:

  Optional annotation column to use as matrix row names.

- normalize_colnames:

  Whether to normalize sample IDs in matrix column names.

## Value

A list with `annotation` and `matrix`.
