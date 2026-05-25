# Fit reference-cohort protein CPSA models

Protein convenience alias for the modality-neutral
[`ageTMP_fit_reference_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_cpsa.md)
engine.

## Usage

``` r
ageTMP_fit_reference_protein_cpsa(
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

A data frame with coefficients, p-values, BY-adjusted FDR values, and
signed combined age-stratum statistics.
