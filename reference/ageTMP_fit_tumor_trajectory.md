# Fit tumor age-dependent molecular trajectories

Fit sex-stratified age-dependent tumor molecular trajectories from a
public feature-by-sample matrix and sample metadata. This is the
tumor-only AD-TMP fitting step used by manuscript heatmap-style
trajectory panels such as Figure 2A.

## Usage

``` r
ageTMP_fit_tumor_trajectory(
  tumor_mat,
  tumor_metadata,
  features = rownames(tumor_mat),
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  center_age_range = c(0, 26),
  fit_age_range = NULL,
  pre_scale = FALSE,
  span = 1.5,
  adaptive_span = FALSE,
  min_span = 0.5,
  max_span = 3,
  span_step = 0.1,
  prediction_ages = NULL,
  prediction_sample_ids = NULL,
  ci_level = 0.95,
  n_cores = 1,
  progress = FALSE
)
```

## Arguments

- tumor_mat:

  Tumor feature-by-sample matrix.

- tumor_metadata:

  Tumor sample metadata.

- features:

  Features to model. Defaults to all matrix rows.

- tumor_sample_col:

  Tumor metadata sample ID column.

- tumor_age_col:

  Tumor metadata age column.

- tumor_sex_col:

  Tumor metadata sex column.

- center_age_range:

  Age range used to center and scale each feature.

- fit_age_range:

  Optional age range used for loess fitting.

- pre_scale:

  Whether to row-center and row-scale the tumor matrix before the
  age-range centering step. This mirrors manuscript trajectory scripts
  that first standardized the complete feature matrix, then standardized
  again within the modeled age range.

- span:

  Loess span. May be a single number or a span data frame accepted by
  the normal/tumor trajectory functions.

- adaptive_span:

  Recompute feature/sex-specific spans by GCV.

- min_span:

  Minimum span for adaptive selection.

- max_span:

  Maximum span for adaptive selection.

- span_step:

  Span grid step for adaptive selection.

- prediction_ages:

  Optional ages where trajectories are predicted. If `NULL`, predictions
  are made at all harmonized tumor sample ages.

- prediction_sample_ids:

  Optional IDs corresponding to `prediction_ages`.

- ci_level:

  Confidence level for fitted trajectories.

- n_cores:

  Number of parallel worker processes for feature/sex fits on Unix-like
  systems. Use `1` for serial execution.

- progress:

  Whether to print simple progress messages.

## Value

A data frame with one row per feature, sex, and prediction age.

## Details

The function intentionally mirrors key details from the manuscript tumor
trajectory workflow. Tumor sample IDs are harmonized with
[`ageTMP_normalize_sample_ids()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_normalize_sample_ids.md),
rows are centered/scaled using the requested `center_age_range`,
optional fitting can be restricted with `fit_age_range`, and
sex-stratified [`stats::loess()`](https://rdrr.io/r/stats/loess.html)
models are fit with base loess defaults.

For exact manuscript reproduction, the original protein workflow used
feature/sex-specific adaptive loess spans selected by generalized
cross-validation over a span grid. Set `adaptive_span = TRUE` to
recompute those spans from the public tumor matrix; this requires the
suggested `locfit` package and can be slow for thousands of proteins. A
single numeric `span`, or a data frame with `feature`, `sex`, `tissue`,
and `span` columns, can also be supplied when spans are known or when a
faster reconstruction is desired.
