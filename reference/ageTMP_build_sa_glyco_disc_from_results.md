# Assemble manuscript STable4 SA-Glyco-Disc columns from glyco CPSA results

Assemble manuscript STable4 SA-Glyco-Disc columns from glyco CPSA
results

## Usage

``` r
ageTMP_build_sa_glyco_disc_from_results(
  female_adj,
  male_adj,
  female_noadj,
  male_noadj,
  membrane_features = NULL
)
```

## Arguments

- female_adj, male_adj, female_noadj, male_noadj:

  Glyco CPSA result data frames, usually the `hope.x00` objects from the
  four manuscript glyco ADO result workspaces.

- membrane_features:

  Optional character vector, data frame, or package membrane-feature
  table identifying membrane-localized glycopeptides.

## Value

A data frame with the published `SA-Glyco-Disc` columns.
