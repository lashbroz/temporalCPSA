# Add age-class columns to a clinical data frame

Add age-class columns to a clinical data frame

## Usage

``` r
ageTMP_add_age_class(
  clinical,
  age_col = "age",
  class_col = "age_class",
  range_col = "age_class_range",
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

- age_col:

  Column containing numeric age.

- class_col:

  Name of the derived age-class column to create.

- range_col:

  Optional name of the derived age-range column to create. Set to `NULL`
  to omit the range column.

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

`clinical` with derived age-class columns appended.
