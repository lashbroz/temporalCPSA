# temporalCPSA Quick Start

`temporalCPSA` supports Cross-Population Survival Analysis (CPSA) and temporal
molecular profiling workflows across independent cohorts.

After installation:

```r
library(temporalCPSA)
help(package = "temporalCPSA")
?temporalCPSA
```

## Core CPSA Workflow

Users provide a feature-by-sample matrix and a clinical/survival table. Age
classes are optional: users may bring their own age-class structure, use
biologically motivated classes, explore age classes from AD-TMP diagnostics, or
omit age classes entirely.

```r
spec <- ageTMP_cpsa_spec(
  time_col = "days",
  event_col = "os.status",
  age_class_col = "age_class",
  age_reference = "ADULT",
  strata = c("PED", "ADO", "YA"),
  base_covariates = c("age", "Grade"),
  adjustment_covariates = character()
)

fit <- ageTMP_fit_reference_cpsa(
  feature_matrix = feature_matrix,
  survival_data = survival_data,
  spec = spec
)
```

For a standard Cox model without age-class interactions:

```r
standard_spec <- ageTMP_cpsa_spec(
  age_class_col = NULL,
  age_reference = NULL,
  strata = NULL,
  include_all_combined_test = FALSE,
  base_covariates = c("age", "Grade"),
  adjustment_covariates = character()
)
```

## Optional Age-Class Exploration

Age-class estimation is an exploratory module, not a decisive automated step in
the CPSA pipeline. It can help assess whether age-dependent temporal molecular
profiles support biologically coherent age intervals.

```r
diagnostic <- ageTMP_segment_diagnostic_matrix(
  tmp_matrix = combined_tmp,
  sample_age = setNames(clinical$age, clinical$id),
  k_values = 2:8,
  min_n = 10
)

ageTMP_plot_segment_diagnostic(
  diagnostic,
  selected_k = 5,
  suggested_cutpoints = c(15, 26, 40, 62),
  suggested_labels = c("PED", "ADO", "YA", "ADULT", "SEN")
)
```

## Useful Help Pages

```r
?ageTMP_cpsa_spec
?ageTMP_fit_reference_cpsa
?ageTMP_predict_tumor_trajectory_matrix
?ageTMP_rank_trajectory_sd
?ageTMP_segment_age_classes
?ageTMP_plot_segment_diagnostic
```

The source README is available at <https://github.com/lashbroz/temporalCPSA>.
