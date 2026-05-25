# Predict tumor age trajectories as a feature matrix

This helper wraps
[`ageTMP_fit_tumor_trajectory()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_tumor_trajectory.md)
for CPSA workflows where the age-dependent tumor trajectory is evaluated
at the ages of an external reference cohort. The result is a numeric
feature-by-sample matrix that can be passed directly to
[`ageTMP_fit_reference_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_cpsa.md).

## Usage

``` r
ageTMP_predict_tumor_trajectory_matrix(
  tumor_mat,
  tumor_metadata,
  prediction_metadata,
  features = rownames(tumor_mat),
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  prediction_sample_col = "id",
  prediction_age_col = "age",
  prediction_sex_col = "Gender",
  center_age_range = c(0, 62),
  fit_age_range = c(0, 80),
  pre_scale = FALSE,
  span = 1.5,
  adaptive_span = FALSE,
  min_span = 0.5,
  max_span = 2,
  span_step = 0.1,
  ci_level = 0.95,
  n_cores = 1,
  progress = FALSE,
  prediction_scope = c("matching_sex", "all_samples"),
  return_trajectory = TRUE
)
```

## Arguments

- tumor_mat:

  Tumor feature-by-sample matrix.

- tumor_metadata:

  Tumor sample metadata.

- prediction_metadata:

  Data frame containing sample IDs, ages, and sex labels for the cohort
  where trajectory values should be predicted.

- features:

  Features to model. Defaults to all matrix rows.

- tumor_sample_col:

  Tumor metadata sample ID column.

- tumor_age_col:

  Tumor metadata age column.

- tumor_sex_col:

  Tumor metadata sex column.

- prediction_sample_col:

  Sample ID column in `prediction_metadata`.

- prediction_age_col:

  Age column in `prediction_metadata`.

- prediction_sex_col:

  Sex column in `prediction_metadata`.

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

- ci_level:

  Confidence level for fitted trajectories.

- n_cores:

  Number of parallel worker processes for feature/sex fits on Unix-like
  systems. Use `1` for serial execution.

- progress:

  Whether to print simple progress messages.

- prediction_scope:

  Whether each sex-specific trajectory should be predicted only for
  samples with the matching sex label (`"matching_sex"`) or for all
  samples in `prediction_metadata` (`"all_samples"`). The latter is
  useful when constructing sex-stratified trajectory rows over a common
  sample grid.

- return_trajectory:

  Whether to include the long-form trajectory table in the result. Set
  to `FALSE` for large CPSA matrices where only the feature-by-sample
  prediction matrix is needed.

## Value

A named list with `matrix` and the long-form `trajectory` table.
