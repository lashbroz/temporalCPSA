# Load reference-cohort mutation covariates

Public-facing alias for
[`ageTMP_load_cdisc_mutation()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_mutation.md).
The returned matrix is suitable for reference-cohort CPSA adjustment
covariates.

## Usage

``` r
ageTMP_load_reference_mutation(
  data_dir = "data",
  genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- genes:

  Character vector of genes to return.

## Value

A numeric sample-by-gene mutation indicator matrix.
