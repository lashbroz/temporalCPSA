# Add manuscript glyco local-FDR columns to a validation result

This is the package analogue of the glyco `vali.x0` post-processing
block in `get_glycosig.R`. It converts a combined p-value column to
one-tailed z scores, fits `locfdr`, and adds `.locfdr.sig`,
`.locfdr.cont`, and `.signed.locfdr.cont` columns for PED and ADO using
the same locfdr threshold.

## Usage

``` r
ageTMP_add_glyco_locfdr_columns(
  fit,
  p_column = "ADO.comb.p",
  classes = c("ADO", "PED"),
  df = 40,
  nulltype = 1,
  pct0 = NULL,
  bre = NULL,
  fdr_cut = 0.1,
  grid = seq(from = 0, to = 1, by = 1e-04),
  main = NULL
)
```

## Arguments

- fit:

  A glyco validation CPSA result data frame.

- p_column:

  Combined p-value column used for the locfdr fit. The manuscript
  no-adjusted glyco blocks use `"ADO.comb.p"` for both PED and ADO.

- classes:

  Combined effect labels to annotate, usually `c("ADO", "PED")`.

- df, nulltype, pct0, bre:

  Arguments passed to
  [`locfdr::locfdr()`](https://rdrr.io/pkg/locfdr/man/locfdr.html).

- fdr_cut:

  Local-FDR cutoff used for the discrete `.locfdr.sig` column.

- grid:

  Numeric grid for the continuous local-FDR lookup.

- main:

  Optional plot title passed to `locfdr`.

## Value

`fit` with glyco locfdr columns added.
