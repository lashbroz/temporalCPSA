# Load a public feature matrix

Load one or more public molecular abundance matrices as
feature-by-sample numeric matrices. Rows are identified by
molecule-level identifiers when available, and duplicate feature rows
can be averaged to produce one row per molecule. For phosphosite data,
uncollapsed matrices use `Site` row IDs, while collapsed matrices use
`ApprovedGeneSymbol` row IDs.

## Usage

``` r
ageTMP_load_feature_matrix(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho"),
  collapse = TRUE,
  row_id = NULL
)
```

## Arguments

- data_dir:

  Path to the public data directory.

- modality:

  Molecular modality or modalities accepted by
  [`ageTMP_load_molecular()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_molecular.md).
  A single modality returns one matrix; multiple modalities return a
  named list of matrices.

- collapse:

  Whether to average rows with the same feature identifier.

- row_id:

  Annotation column to use as the feature identifier. If `NULL`, a
  modality-aware default is used.

## Value

A numeric feature-by-sample matrix, or a named list of matrices when
multiple modalities are requested.
