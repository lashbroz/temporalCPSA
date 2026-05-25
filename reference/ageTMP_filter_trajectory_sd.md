# Filter trajectory features by standard deviation

Apply an explicit dynamic-range filter to trajectory SD rankings. The
default keeps features with `sd >= 0.15`, matching the omit-list
threshold used in the manuscript `rank_df_rev.RData` workflow.

## Usage

``` r
ageTMP_filter_trajectory_sd(
  sd_rank,
  sd_min = 0.15,
  group_cols = c("sex"),
  keep = c("all_groups", "any_group", "per_group")
)
```

## Arguments

- sd_rank:

  Data frame returned by
  [`ageTMP_rank_trajectory_sd`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_rank_trajectory_sd.md).

- sd_min:

  Minimum standard deviation required to keep a feature.

- group_cols:

  Columns defining groups that must pass the threshold.

- keep:

  One of `"all_groups"`, `"any_group"`, or `"per_group"`.

## Value

A character vector of kept features for `"all_groups"` or `"any_group"`,
and a filtered data frame for `"per_group"`.
