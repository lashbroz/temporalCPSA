# Add glyco local-FDR columns with manuscript sex-specific settings

Add glyco local-FDR columns with manuscript sex-specific settings

## Usage

``` r
ageTMP_add_glyco_locfdr_by_manuscript_settings(
  fit,
  sex = c("Female", "Male"),
  adjusted = FALSE,
  ...
)
```

## Arguments

- fit:

  A glyco validation CPSA result data frame.

- sex:

  One of `"Female"` or `"Male"`.

- adjusted:

  Whether the glyco result is protein-adjusted.

- ...:

  Additional arguments passed to
  [`ageTMP_add_glyco_locfdr_columns()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_glyco_locfdr_columns.md).

## Value

`fit` with glyco locfdr columns added.
