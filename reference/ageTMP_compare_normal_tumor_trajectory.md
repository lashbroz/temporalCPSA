# Compare normal and tumor age trajectories

This function implements the Figure 2-style normal/tumor trajectory
comparison using public tumor molecular data and package-stored normal
reference data. Normal values are first adjusted for technical
covariates such as pH, PMI, and ethnicity, then normal and tumor age
trajectories are fit separately by sex.

## Usage

``` r
ageTMP_compare_normal_tumor_trajectory(
  tumor_mat,
  tumor_metadata,
  normal_mat,
  normal_metadata,
  features,
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  normal_sample_col = "ID",
  normal_age_col = "Age",
  normal_sex_col = "Gender",
  normal_covariates = c("pH", "PMI", "Ethnicity"),
  center_age_range = c(0, 26),
  span = 1.5,
  adaptive_span = FALSE,
  tumor_min_span = 0.5,
  tumor_max_span = 3,
  normal_min_span = 1,
  normal_max_span = 3,
  span_step = 0.1,
  prediction_ages = NULL,
  prediction_sample_ids = NULL,
  ci_level = 0.95
)
```

## Arguments

- tumor_mat:

  Tumor feature-by-sample matrix.

- tumor_metadata:

  Tumor sample metadata.

- normal_mat:

  Normal/reference feature-by-sample matrix.

- normal_metadata:

  Normal/reference sample metadata.

- features:

  Features to plot/model.

- tumor_sample_col:

  Tumor metadata sample ID column.

- tumor_age_col:

  Tumor metadata age column.

- tumor_sex_col:

  Tumor metadata sex column.

- normal_sample_col:

  Normal metadata sample ID column.

- normal_age_col:

  Normal metadata age column.

- normal_sex_col:

  Normal metadata sex column.

- normal_covariates:

  Covariates used to adjust normal data.

- center_age_range:

  Age range used to center/scale tumor and normal data.

- span:

  Loess span for trajectory fitting. This can be a single numeric value
  used for all fits, a list with `Normal` and `Tumor` entries, or a data
  frame with columns `feature`, `sex`, `tissue`, and `span`.

- adaptive_span:

  Recompute feature/sex/tissue-specific spans by GCV. This mirrors the
  adaptive protein trajectory run used for Figure 2F/G.

- tumor_min_span:

  Minimum span for adaptive tumor trajectories.

- tumor_max_span:

  Maximum span for adaptive tumor trajectories.

- normal_min_span:

  Minimum span for adaptive normal trajectories.

- normal_max_span:

  Maximum span for adaptive normal trajectories.

- span_step:

  Span grid step for adaptive selection.

- prediction_ages:

  Optional numeric vector of ages where trajectories should be
  predicted. If `NULL`, trajectories are predicted at the tumor sample
  ages for the sex being modeled. For manuscript Figure 2C reproduction,
  pass all tumor sample ages `<= 50` so each sex-stratified curve is
  evaluated on the same age support used by the original `tn.df` object.

- prediction_sample_ids:

  Optional sample IDs corresponding to `prediction_ages`. Supplying
  these preserves manuscript heatmap columns when multiple samples share
  the same age.

- ci_level:

  Confidence level for fitted trajectories.

## Value

A data frame with fitted values, standard errors, and confidence
intervals for normal and tumor trajectories at tumor sample ages or at
`prediction_ages` when provided.

## Details

Several defaults and preprocessing steps intentionally preserve the
original manuscript trajectory workflow rather than imposing a more
general modeling API.

- Loess fits use base
  [`stats::loess()`](https://rdrr.io/r/stats/loess.html) defaults,
  matching the original trajectory scripts. In particular, the package
  does not use `loess.control(surface = "direct")` for manuscript
  reproduction.

- Normal/reference values are adjusted with
  `score ~ pH + PMI + Ethnicity` when those covariates are supplied.
  Missing normal pH is imputed before sex-stratified fitting; ethnicity
  label `H` is combined with `C`; and `C` is used as the reference
  level.

- Tumor matrices should be centered/scaled using the manuscript
  centering range, e.g. `center_age_range = c(0, 50)` for the Figure 2C
  protein trajectories.

- The published Figure 2C protein workflow fit tumor curves using
  samples with `age <= 80` (`age.cut2 = 80` in `proteo_tadj50.R`). This
  filtering is performed by the reproduction script before calling this
  function.

- For Figure 2C reproduction, each sex-stratified model is predicted on
  the common set of all tumor sample ages `<= 50`, not only on ages from
  the sex being modeled. Supply those ages through `prediction_ages`.

- Manuscript panels used feature/sex/tissue-specific adaptive loess
  spans. Supply them as a data frame with columns `feature`, `sex`,
  `tissue`, and `span`, or set `adaptive_span = TRUE` to recompute them
  from the supplied public source data.

- Confidence intervals are returned as `fit +/- qnorm(0.975) * se` when
  `ci_level = 0.95`.

These details were validated against the original `tn_df.RData` and
`protein_tadj50_list.RData` objects for the Figure 2C proteins `CNTN1`,
`MAPT`, and `L1CAM`; with the manuscript settings, generated `fit`,
`ci_lower`, and `ci_upper` values match the legacy `tn.df$value`,
`tn.df$low`, and `tn.df$hi` to floating-point tolerance.
