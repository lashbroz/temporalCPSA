# Load reference-cohort protein matrix

Public-facing alias for
[`ageTMP_load_cdisc_protein_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_protein_matrix.md).
The manuscript implementation currently supports protein CPSA first; the
naming leaves room for RNA, glyco, phosphosite, and other feature
matrices to use the same reference-cohort vocabulary.

## Usage

``` r
ageTMP_load_reference_protein_matrix(data_dir = "data", collapse = TRUE)
```

## Arguments

- data_dir:

  Path to the public data directory.

- collapse:

  Whether to average rows with the same gene symbol.

## Value

A numeric gene-by-sample protein matrix.
