# Build the supplementary age-class table

Build the `cDisc_AgeClass`-style supplementary table from clinical ages
and, when available, the TMP matrix used for silhouette scoring.
Age-class labels and ranges are derived from numeric ages. `sil_width`
is computed only when a feature-by-sample TMP matrix is supplied.

## Usage

``` r
ageTMP_build_age_class_supplement(
  clinical,
  tmp_matrix = NULL,
  sample_ids = NULL,
  id_col = "id",
  age_col = "age",
  breaks = c(0, 15, 26, 40, 62, 80),
  labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
)
```

## Arguments

- clinical:

  Clinical data frame.

- tmp_matrix:

  Optional feature-by-sample TMP matrix used to compute sample
  silhouette widths.

- sample_ids:

  Optional sample IDs to include. If `NULL`, all clinical samples with
  non-missing age are used, or the intersection with
  `colnames(tmp_matrix)` when a TMP matrix is supplied.

- id_col:

  Clinical sample ID column.

- age_col:

  Clinical age column.

- breaks:

  Numeric breakpoints passed to
  [`base::cut()`](https://rdrr.io/r/base/cut.html).

- labels:

  Age-class labels. Length must be `length(breaks) - 1`.

- right:

  Whether intervals are closed on the right.

- include_lowest:

  Whether the lowest age should be included in the first interval.

- ordered:

  Whether to return an ordered factor.

## Value

A data frame with `id`, `age.class`, `age.class.range`, and `sil_width`.
