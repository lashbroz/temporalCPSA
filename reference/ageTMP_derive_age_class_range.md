# Derive age-class interval labels from age

Return the interval labels corresponding to
[`ageTMP_derive_age_class()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_derive_age_class.md).
This is useful when reproducing manuscript tables that show both a
semantic class label, such as `PED`, and the numeric range, such as
`[0,15]`.

## Usage

``` r
ageTMP_derive_age_class_range(
  age,
  breaks = c(0, 15, 26, 40, 62, 80),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
)
```

## Arguments

- age:

  Numeric age vector.

- breaks:

  Numeric breakpoints passed to
  [`base::cut()`](https://rdrr.io/r/base/cut.html).

- right:

  Whether intervals are closed on the right.

- include_lowest:

  Whether the lowest age should be included in the first interval.

- ordered:

  Whether to return an ordered factor.

## Value

A factor of interval labels.
