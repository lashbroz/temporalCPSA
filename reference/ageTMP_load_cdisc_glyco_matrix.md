# Load public cDisc glycopeptide data

Unlike gene-level protein/RNA loaders, glycopeptide survival analyses
operate at the glycopeptide row level. This loader preserves
`Gene.Sequence` row identifiers and returns the full public glycopeptide
annotation.

## Usage

``` r
ageTMP_load_cdisc_glyco_matrix(data_dir = "data")
```

## Arguments

- data_dir:

  Path to the public data directory.

## Value

A list with `matrix` and `annotation` elements.
