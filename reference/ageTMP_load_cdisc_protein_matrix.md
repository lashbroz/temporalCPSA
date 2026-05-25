# Load manuscript reference-cohort protein data as a gene-by-sample matrix

Load manuscript reference-cohort protein data as a gene-by-sample matrix

## Usage

``` r
ageTMP_load_cdisc_protein_matrix(data_dir = "data", collapse = TRUE)
```

## Arguments

- data_dir:

  Path to the public data directory.

- collapse:

  Whether to average rows with the same gene symbol.

## Value

A numeric gene-by-sample matrix.
