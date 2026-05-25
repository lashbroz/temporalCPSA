# Build source-derived protein-adjusted glycopeptide matrix

This reproduces the manuscript adjusted glycopeptide matrix by
residualizing each glycopeptide abundance vector against the matched
gene-level protein abundance vector over shared cDisc samples.

## Usage

``` r
ageTMP_build_cdisc_glyco_adjusted_matrix(
  data_dir = "data",
  glyco = NULL,
  protein_matrix = NULL
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- glyco:

  Optional output from
  [`ageTMP_load_cdisc_glyco_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_glyco_matrix.md).

- protein_matrix:

  Optional gene-by-sample protein matrix.

## Value

A list with `matrix` and `annotation` elements.
