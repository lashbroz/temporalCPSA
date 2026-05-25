# Load manuscript trajectory QC omit lists

The manuscript CPSA workflow uses sex-specific omit lists for features
with very flat age-dependent trajectory signal. These omit lists are
derived from fitted trajectory matrices by ranking the row-wise standard
deviation and using `sd < 0.15` as the default low-dynamic-range
threshold.

## Usage

``` r
ageTMP_load_trajectory_qc(rdata_path, type = "protein")
```

## Arguments

- rdata_path:

  Path to a `rank_df_rev.RData`-style object.

- type:

  Molecular data type to keep, such as `"protein"`.

## Value

A list with the rank data frame and male/female omit vectors.
