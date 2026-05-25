# Load manuscript reference-cohort mutation calls for CPSA modeling

The manuscript protein CPSA workflow uses the reference-cohort mutation
modeling matrix `cDisc_mutation_10192023.tsv`, which includes both local
discovery samples and external reference-labeled samples present in the
survival cohort. If that file is unavailable, the function falls back to
the `Disc_Mutation` sheet in `STable1.xlsx`, but that reduced
supplementary view does not contain the full mutation covariate set used
by the manuscript CPSA models.

## Usage

``` r
ageTMP_load_cdisc_mutation(
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

A numeric sample-by-gene matrix with mutation indicators.
