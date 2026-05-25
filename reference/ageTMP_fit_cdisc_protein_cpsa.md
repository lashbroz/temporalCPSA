# Fit reference-cohort protein CPSA Cox models

This fits the reference-cohort CPSA model used for manuscript protein
features: a Cox model with scaled protein abundance,
protein-by-age-class interactions, age, tumor covariates, and mutation
covariates. Combined age-stratum p-values are obtained by dropping the
stratum-specific protein term from the full model.

## Usage

``` r
ageTMP_fit_cdisc_protein_cpsa(
  protein_matrix,
  survival_data,
  genes = rownames(protein_matrix),
  age_reference = "ADULT",
  mutation_genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
)
```

## Arguments

- protein_matrix:

  Numeric gene-by-sample protein matrix.

- survival_data:

  Prepared survival covariates from
  [`ageTMP_prepare_reference_cpsa_survival()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_prepare_reference_cpsa_survival.md).

- genes:

  Optional subset of genes to fit.

- age_reference:

  Reference age class for the interaction model.

- mutation_genes:

  Mutation gene names without `_mut`.

## Value

A data frame with reference-cohort CPSA statistics.
