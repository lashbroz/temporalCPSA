# Normalize HOPE AYA sample identifiers

This mirrors the repeated `get_id()` helper used in the paper scripts.

## Usage

``` r
ageTMP_normalize_sample_ids(x)
```

## Arguments

- x:

  Character vector of sample identifiers.

## Value

A character vector with common prefixes removed and periods converted to
dashes.
