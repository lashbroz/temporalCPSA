# Derive age-class labels from age

Convert numeric ages to ordered age-class labels using explicit
breakpoints. The default breakpoints reproduce the manuscript five-class
age grouping: PED, ADO, YA, ADULT, and SEN.

## Usage

``` r
ageTMP_derive_age_class(
  age,
  breaks = c(0, 15, 26, 40, 62, 80),
  labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
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

- labels:

  Age-class labels. Length must be `length(breaks) - 1`.

- right:

  Whether intervals are closed on the right.

- include_lowest:

  Whether the lowest age should be included in the first interval.

- ordered:

  Whether to return an ordered factor.

## Value

An ordered factor of age-class labels.
