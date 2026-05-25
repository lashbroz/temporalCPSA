# Load a public molecular data table

Load a public molecular data table

## Usage

``` r
ageTMP_load_molecular(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho", "mutation", "full_mutation")
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- modality:

  One of `"protein"`, `"rna"`, `"glyco"`, `"phospho"`, `"mutation"`, or
  `"full_mutation"`.

## Value

A data frame containing feature annotation columns followed by sample
columns.
