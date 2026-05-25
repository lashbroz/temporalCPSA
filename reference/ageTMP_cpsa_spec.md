# Create a CPSA model specification

CPSA is feature-matrix agnostic: each row of `feature_matrix` is a
molecular feature or score, and each column is a sample. This
specification object makes the survival outcome, age-stratum
interaction, and clinical/molecular adjustment model explicit so users
can adapt CPSA to different applications.

## Usage

``` r
ageTMP_cpsa_spec(
  time_col = "days",
  event_col = "os.status",
  feature_col = "feature_score",
  age_class_col = "age_class_new",
  age_reference = "ADULT",
  strata = c("PED", "ADO", "YA"),
  include_all_combined_test = TRUE,
  base_covariates = c("age", "Grade", "Cortical", "Midline"),
  adjustment_covariates = c("IDH1_mut", "TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut"),
  fdr_method = "BY",
  remove_sparse_covariates = TRUE,
  sparse_covariates = NULL,
  sparse_positive_level = 1,
  sparse_min_count = 3,
  sparse_count_method = c("numeric", "factor"),
  scale_features = TRUE,
  engine = c("coxph", "coxph.fit")
)
```

## Arguments

- time_col:

  Column in `survival_data` containing survival time.

- event_col:

  Column in `survival_data` containing event status.

- feature_col:

  Temporary column name used for the current feature.

- age_class_col:

  Column containing age-class labels. Set to `NULL` with
  `age_reference = NULL` and `strata = NULL` to fit a standard feature
  Cox model without feature-by-age-class interactions.

- age_reference:

  Reference level for `age_class_col`.

- strata:

  Age strata to test by dropping feature terms.

- include_all_combined_test:

  Whether to add an omnibus `ALL.comb.*` likelihood-ratio test that
  drops the main feature and all feature-by-strata interaction terms.
  Legacy two-level glyco ADO scripts did not create this column.

- base_covariates:

  Clinical covariates included in every model.

- adjustment_covariates:

  Additional covariates, such as mutation calls.

- fdr_method:

  Multiple-testing adjustment method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- remove_sparse_covariates:

  Whether to drop binary covariates with very few positive values in the
  current model frame.

- sparse_covariates:

  Optional character vector limiting which covariates are eligible for
  sparse binary checks. If `NULL`, only covariates with no more than two
  observed values are checked.

- sparse_positive_level:

  Positive value used for sparse binary checks.

- sparse_min_count:

  Drop sparse covariates with counts less than or equal to this value.

- sparse_count_method:

  How to count sparse positive values. `"numeric"` coerces
  binary/logical covariates to numeric before counting. `"factor"`
  preserves the legacy manuscript behavior of counting
  `table(factor(x, levels = c(0, 1)))[2]`, which drops logical
  covariates such as `Cortical` and `Midline` in the glyco ADO scripts.

- scale_features:

  Whether to row-center and row-scale features inside the model engine.
  Set to `FALSE` when the supplied matrix has already been scaled in a
  manuscript-specific way.

- engine:

  Cox model fitting engine. `"coxph"` uses the formula-based
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  path and is the default for manuscript replication. `"coxph.fit"` uses
  the same model matrix with
  [`survival::coxph.fit()`](https://rdrr.io/pkg/survival/man/agreg.fit.html)
  to reduce formula overhead in larger feature scans.

## Value

A list describing the CPSA Cox model.
