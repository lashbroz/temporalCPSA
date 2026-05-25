# Return manuscript glyco local-FDR settings

Return manuscript glyco local-FDR settings

## Usage

``` r
ageTMP_glyco_locfdr_settings(sex = c("Female", "Male"), adjusted = FALSE)
```

## Arguments

- sex:

  One of `"Female"` or `"Male"`.

- adjusted:

  Whether the glyco result is protein-adjusted.

## Value

A list of arguments suitable for
[`ageTMP_add_glyco_locfdr_columns()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_glyco_locfdr_columns.md).
