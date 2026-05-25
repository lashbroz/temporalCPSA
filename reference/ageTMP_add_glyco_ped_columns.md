# Add manuscript PED columns to a glyco ADO CPSA result

The glyco supplementary table was generated from PED/ADO-only Cox
models. The likelihood-ratio p-value is therefore shared by PED and ADO,
while the PED sign uses the main feature coefficient and the ADO sign
uses the feature-plus-interaction coefficient already present in
`ADO.comb.*`.

## Usage

``` r
ageTMP_add_glyco_ped_columns(fit)
```

## Arguments

- fit:

  A glyco CPSA result data frame with `ADO.comb.p`, `ADO.comb.fdr`, and
  `score2.coef`.

## Value

`fit` with `PED.comb.*` columns added.
