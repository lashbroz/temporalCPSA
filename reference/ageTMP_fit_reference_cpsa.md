# Fit reference-cohort CPSA Cox models

This is the modality-neutral CPSA model engine. It accepts any numeric
feature-by-sample matrix, including protein abundance, RNA expression,
glycopeptide abundance, phosphosite abundance, pathway scores, or other
molecular feature scores.

## Usage

``` r
ageTMP_fit_reference_cpsa(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  spec = ageTMP_cpsa_spec(),
  trajectory_sd = NULL,
  trajectory_sd_min = NULL,
  trajectory_sd_feature_col = "feature",
  trajectory_sd_value_col = "sd",
  trajectory_sd_group_cols = NULL,
  trajectory_sd_keep = c("any_group", "all_groups"),
  n_cores = 1,
  progress = FALSE,
  progress_every = 500
)
```

## Arguments

- feature_matrix:

  Numeric feature-by-sample matrix.

- survival_data:

  Prepared survival covariates.

- features:

  Optional subset of features to fit.

- spec:

  CPSA model specification from
  [`ageTMP_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_cpsa_spec.md).

- trajectory_sd:

  Optional trajectory dynamic-range object used to screen features
  before CPSA fitting. Accepted inputs are a named numeric vector of
  per-feature trajectory standard deviations, a feature-by-age/sample
  matrix whose row-wise standard deviation should be used, or a data
  frame such as the output of
  [`ageTMP_rank_trajectory_sd()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_rank_trajectory_sd.md).

- trajectory_sd_min:

  Minimum trajectory standard deviation required to keep a feature. Use
  `NULL` to disable trajectory-SD screening.

- trajectory_sd_feature_col:

  Feature column in a `trajectory_sd` data frame.

- trajectory_sd_value_col:

  Standard-deviation column in a `trajectory_sd` data frame.

- trajectory_sd_group_cols:

  Optional grouping columns in a `trajectory_sd` data frame. When
  supplied, `trajectory_sd_keep` controls whether a feature must pass in
  all groups or any group.

- trajectory_sd_keep:

  Whether grouped trajectory-SD screening keeps features passing in
  `"any_group"` or `"all_groups"`.

- n_cores:

  Number of parallel worker processes on non-Windows systems.

- progress:

  Whether to print fitting progress.

- progress_every:

  Number of features per progress chunk.

## Value

A data frame with coefficients, p-values, adjusted FDR values, and
signed age-stratum statistics.
