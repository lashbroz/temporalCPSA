# Rank trajectory features by standard deviation

Compute within-group standard deviation of fitted trajectory values.
This makes explicit the manuscript QC idea behind `rank_df_rev.RData`,
where some downstream analyses exclude features with very flat
age-dependent tumor trajectories, commonly using `sd < 0.15` as the omit
threshold.

## Usage

``` r
ageTMP_rank_trajectory_sd(
  trajectory,
  feature_col = "feature",
  value_col = "fit",
  group_cols = c("sex"),
  tissue_col = "tissue",
  tissue = NULL
)
```

## Arguments

- trajectory:

  Long-format trajectory data frame.

- feature_col:

  Column containing feature names.

- value_col:

  Column containing trajectory values, usually `fit`.

- group_cols:

  Columns defining independent rankings.

- tissue_col:

  Optional tissue column used with `tissue`.

- tissue:

  Optional tissue value to retain before ranking, such as `"Tumor"`.

## Value

A data frame with grouping columns, `feature`, `sd`, and `rank_sd`.
