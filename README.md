# temporalCPSA

`temporalCPSA` is an R package implementing the Cross-Population Survival Analysis (CPSA) framework described in Tignor et al., *Proteogenomic analysis of pediatric and AYA high-grade glioma*, enabling reproducible study execution together with generalized workflows for temporal molecular pattern (TMP) analysis, age-dependent molecular trajectory modeling, and cross-population survival analysis across complementary patient cohorts.

## Installation

Install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("lashbroz/temporalCPSA")
```

## What Is CPSA?

Cross-Population Survival Analysis (CPSA) is a framework for assessing
prognostic molecular trajectories across complementary patient cohorts.

CPSA estimates temporal molecular patterns (TMPs) in deeply profiled discovery
cohorts, evaluates their clinical relevance in external clinically annotated
populations, and validates these prognostic patterns back in the discovery
cohort. By leveraging age-associated molecular structure rather than survival
associations observed within a single cohort alone, CPSA supports
cross-population prioritization of prognostic signals that may reflect
conserved biological dynamics across heterogeneous datasets.

## What Can This Package Be Used For?

- Reproducing analyses from Tignor et al., *Proteogenomic analysis of
  pediatric and AYA high-grade glioma*, using the companion
  [Wang Lab manuscript repository](https://github.com/WangLab-MSSM/Pediatric-AYA-high-grade-glioma).
- Using deeply profiled CPTAC or related multi-omic cohorts to inform survival
  modeling in external clinically annotated populations with sparse or absent
  molecular data.
- Leveraging temporal molecular structure to identify prognostic biomarkers
  whose outcome associations reflect conserved biological dynamics across
  complementary patient cohorts.
- Modeling age-associated molecular trajectories to facilitate comparison of age-related
  molecular dynamics within and between cohorts, including tumor versus normal
  and male versus female tissue comparisons.
- Exploration of age structure within cohorts based on clustering of molecular
  abundance profiles across diverse multi-omic data types.

## General CPSA Workflow

The core workflow starts with a feature-by-sample matrix and a reference-cohort
clinical table containing survival time, event status, age class, and
adjustment covariates.

Users may supply age classes, model age continuously, explore age classes from
AD-TMP diagnostics, or omit age classes entirely. The CPSA model specification
controls which choice is used.

Default arguments in several trajectory and CPSA helpers reflect the manuscript
analysis settings; users should review and modify these settings for other
cohorts, molecular platforms, or modeling designs.

## Estimating Temporal Molecular Trajectories

The basic trajectory-estimation helper is
`ageTMP_predict_tumor_trajectory_matrix()`. It fits age-associated molecular
trajectories in a deeply profiled cohort and predicts those trajectory values
at the ages of a target cohort. The target cohort can be the same cohort, an
external clinically annotated cohort, or an age grid used for visualization.

```r
library(temporalCPSA)

tumor_matrix <- your_feature_by_sample_matrix
discovery_clinical <- your_discovery_clinical_table
target_clinical <- your_target_clinical_table

trajectory_fit <- ageTMP_predict_tumor_trajectory_matrix(
  tumor_mat = tumor_matrix,
  tumor_metadata = discovery_clinical,
  prediction_metadata = target_clinical,
  tumor_sample_col = "sample_id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  prediction_sample_col = "sample_id",
  prediction_age_col = "age",
  prediction_sex_col = "sex",
  span = 0.75,
  prediction_scope = "all_samples",
  return_trajectory = TRUE
)

trajectory_matrix <- trajectory_fit$matrix
trajectory_long <- trajectory_fit$trajectory
```

For tumor-versus-normal/reference trajectory comparisons, use
`ageTMP_compare_normal_tumor_trajectory()` when both a tumor cohort and a
normal/reference cohort are available:

```r
trajectory_comparison <- ageTMP_compare_normal_tumor_trajectory(
  tumor_mat = tumor_matrix,
  normal_mat = normal_reference_matrix,
  tumor_metadata = discovery_clinical,
  normal_metadata = normal_reference_clinical,
  tumor_sample_col = "sample_id",
  normal_sample_col = "sample_id",
  tumor_age_col = "age",
  normal_age_col = "age",
  tumor_sex_col = "sex",
  normal_sex_col = "sex",
  features = rownames(tumor_matrix)[1:25],
  span = 0.75
)
```

```r
feature_matrix <- your_feature_matrix
survival_data <- your_survival_data

spec <- ageTMP_cpsa_spec(
  time_col = "days",
  event_col = "os.status",
  age_class_col = "age_class",
  age_reference = "ADULT",
  strata = c("PED", "ADO", "YA"),
  base_covariates = c("age", "Grade", "Cortical", "Midline"),
  adjustment_covariates = c("TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut")
)

fit <- ageTMP_fit_reference_cpsa(
  feature_matrix = feature_matrix,
  survival_data = survival_data,
  spec = spec
)
```

For a standard Cox model without age-class interaction terms, set
`age_class_col`, `age_reference`, and `strata` to `NULL`:

```r
standard_spec <- ageTMP_cpsa_spec(
  age_class_col = NULL,
  age_reference = NULL,
  strata = NULL,
  include_all_combined_test = FALSE,
  base_covariates = c("age", "Grade"),
  adjustment_covariates = character()
)

standard_fit <- ageTMP_fit_reference_cpsa(
  feature_matrix,
  survival_data,
  spec = standard_spec
)
```

When CPSA is applied to age-dependent trajectory features, users should screen
out features with minimal trajectory dynamic range before interpreting survival
associations. The package can apply this screen inside the CPSA fit using a
named trajectory-SD vector, a trajectory prediction matrix, or the output of
`ageTMP_rank_trajectory_sd()`:

```r
trajectory_sd <- ageTMP_rank_trajectory_sd(
  tumor_trajectory,
  feature_col = "feature",
  value_col = "fit",
  group_cols = "trajectory_stratum"
)

fit <- ageTMP_fit_reference_cpsa(
  feature_matrix,
  survival_data,
  spec = spec,
  trajectory_sd = trajectory_sd,
  trajectory_sd_min = 0.15,
  trajectory_sd_group_cols = "trajectory_stratum",
  trajectory_sd_keep = "any_group"
)
```

Here, `trajectory_sd_group_cols` identifies the column or columns defining
trajectory strata for dynamic-range screening. In the manuscript analyses this
was often sex, but the filter itself is generic. With
`trajectory_sd_keep = "any_group"`, a feature is retained if it is dynamic in at
least one trajectory stratum; use stricter settings when a feature should be
dynamic across all modeled strata.

The CPSA model engine computes raw p-values and, by default, Benjamini-Yekutieli
FDR columns with `stats::p.adjust(..., method = "BY")`. Because CPSA is often
applied across correlated molecular features and trajectory-derived scores,
users should inspect empirical test-statistic distributions and consider
whether additional calibration is appropriate for their modeling design.

In the manuscript reference-cohort analyses, protein and RNA significance calls
used a separate local-FDR procedure on transformed p-values. That procedure is
documented in the manuscript figure-generation scripts and is not part of the
general package API.

## Optional Age-Class Exploration

Age-class estimation is optional. Users may supply age classes, model age
continuously, or use `temporalCPSA` diagnostics to explore whether AD-TMP
structure supports age-contiguous intervals.

In the manuscript, age classes were chosen by considering both temporal
molecular clustering and conventional age distributions. In `temporalCPSA`, we
implement a data-driven diagnostic approach: candidate classes are constrained
to ordered age intervals, while the AD-TMP matrix determines which intervals are
best supported.

These diagnostics are intended for cohort exploration and interpretation, not
as a required or decisive step in the CPSA pipeline.

```r
library(temporalCPSA)

clinical <- ageTMP_load_discovery_clinical("data")
feature_matrices <- ageTMP_load_feature_matrix(
  "data",
  modality = c("protein", "rna", "phospho")
)

tmp_matrices <- lapply(feature_matrices, function(feature_matrix) {
  ageTMP_predict_tumor_trajectory_matrix(
    tumor_mat = feature_matrix,
    tumor_metadata = clinical,
    prediction_metadata = clinical,
    tumor_sample_col = "id",
    tumor_age_col = "cDisc_age",
    tumor_sex_col = "cDisc_Gender",
    prediction_sample_col = "id",
    prediction_age_col = "cDisc_age",
    prediction_sex_col = "cDisc_Gender",
    return_trajectory = FALSE
  )$matrix
})

combined_tmp <- ageTMP_combine_tmp_matrices(tmp_matrices)
age_clusters <- ageTMP_cluster_tmp_samples(combined_tmp)
```

To derive exploratory age classes with guaranteed contiguity:

```r
fit <- ageTMP_segment_age_classes(
  tmp_matrix = combined_tmp,
  sample_age = setNames(clinical$cDisc_age, clinical$id),
  k = 5,
  min_n = 10,
  labels = c("PED", "ADO", "YA", "ADULT", "SEN")
)

fit$cutpoints
fit$segment_summary
```

To review the ordered segmentation landscape, draw the green diagnostic heatmap
and overlay optimizer-selected, manuscript-style, or investigator-selected
cutpoints.

```r
diagnostic <- ageTMP_segment_diagnostic_matrix(
  tmp_matrix = combined_tmp,
  sample_age = setNames(clinical$cDisc_age, clinical$id),
  k_values = 2:8,
  min_n = 10
)

ageTMP_plot_segment_diagnostic(
  diagnostic,
  selected_k = 5,
  class_labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
  suggested_cutpoints = c(15, 26, 40, 62),
  suggested_labels = c("PED", "ADO", "YA", "ADULT", "SEN")
)
```

Simple helpers such as `ageTMP_derive_age_class()` remain available for applying
predefined cutpoints.

## Manuscript Reproduction

The package also supports manuscript reproduction from public source files:

Paper figure-generation scripts live in the companion manuscript repository:
[`WangLab-MSSM/Pediatric-AYA-high-grade-glioma`](https://github.com/WangLab-MSSM/Pediatric-AYA-high-grade-glioma).

```r
# Example manuscript workflow pieces
temporalCPSA::ageTMP_data_sources("data")
temporalCPSA::ageTMP_load_normal_reference()
```

Preferred public-facing functions use reference-cohort language:

```r
clinical <- ageTMP_load_reference_clinical("data")
mutation <- ageTMP_load_reference_mutation("data")
feature_matrix <- ageTMP_load_reference_protein_matrix("data")

survival_data <- ageTMP_prepare_reference_cpsa_survival(
  clinical = clinical,
  mutation = mutation,
  protein_sample_ids = colnames(feature_matrix),
  sex = "Female",
  mutation_na = "preserve"
)

manuscript_spec <- ageTMP_cpsa_spec(
  base_covariates = c("age", "Grade", "Cortical", "Midline"),
  adjustment_covariates = c("TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut")
)

manuscript_fit <- ageTMP_fit_reference_cpsa(
  feature_matrix,
  survival_data,
  spec = manuscript_spec
)
```

Manuscript-specific wrappers such as `ageTMP_build_sa_protein_cdisc()` are kept
for direct table reproduction because the published table columns use the
`cdisc` label.

## Terminology

The package uses general cohort language in its public API:

- **discovery cohort**: the smaller, deeply profiled tumor cohort used to define
  age-dependent tumor molecular phenotypes and candidate molecular programs;
- **reference cohort**: an external or larger cohort with compatible molecular
  features, survival outcomes, clinical covariates, and optional mutation or
  other adjustment covariates;
- **feature matrix**: a numeric feature-by-sample matrix. Rows may represent
  protein abundance, RNA expression, glycopeptide abundance, phosphosite
  abundance, methylation features, pathway scores, or other molecular scores;
- **CPSA**: Cross-Population Survival Analysis, the downstream survival modeling
  framework used to evaluate whether discovery-cohort molecular programs show
  age- or sex-stratified survival associations in a reference cohort.

The manuscript source files and published supplementary tables retain labels
such as `cDisc` in file names and column names. In package documentation, `cDisc`
should be read as the manuscript's reference cohort label, not as a requirement
for general use.

Normal/reference DLPFC developmental data used for trajectory analyses are
documented package extdata derived from the DEveLopmental Trajectory Atlas
(DELTA), PMID: 30518843, available at <http://amp.pharm.mssm.edu/DELTA>.

## Reproducibility Principle

The complete manuscript figure-generation scripts are maintained in the
Wang Lab reproducibility repository:
[`WangLab-MSSM/Pediatric-AYA-high-grade-glioma`](https://github.com/WangLab-MSSM/Pediatric-AYA-high-grade-glioma).

Reproducible workflows should read directly from documented source files
whenever possible and avoid hidden `.RData` objects or manually generated
intermediate TSV files. Serialized package data are reserved for documented
external reference data that are part of the analysis provenance.
