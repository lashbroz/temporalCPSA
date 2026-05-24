# temporalCPSA

`temporalCPSA` is an R package implementing the Cross-Population Survival Analysis (CPSA) framework described in Tignor et al., *Proteogenomic analysis of pediatric and AYA high-grade glioma*, supporting both reproducible execution of published study analyses and generalized workflows for integrating temporal molecular profiles with clinical outcome data across heterogeneous multi-omic cohorts.

## Installation

```r
install.packages("remotes")
remotes::install_github("lashbroz/temporalCPSA")
```

## General CPSA Workflow

The core workflow starts with a feature-by-sample matrix and a reference-cohort
clinical table. Rows of the feature matrix are molecular features or trajectory
scores, and columns are sample IDs. The clinical table supplies survival time,
event status, age class, and adjustment covariates.

Users may bring their own age-class structure, use age classes defined by a
prior biological or clinical rationale, explore age classes from AD-TMP
diagnostics, or omit age classes entirely. The CPSA model specification controls
which of these choices is used for a given analysis.

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
  group_cols = "sex"
)

fit <- ageTMP_fit_reference_cpsa(
  feature_matrix,
  survival_data,
  spec = spec,
  trajectory_sd = trajectory_sd,
  trajectory_sd_min = 0.15,
  trajectory_sd_group_cols = "sex",
  trajectory_sd_keep = "any_group"
)
```

## Optional Age-Class Exploration

Age-class estimation is optional. Users may supply their own age classes, model
age continuously, or use `temporalCPSA` diagnostics to explore whether AD-TMP
structure supports age-contiguous intervals.

In the manuscript, age classes were chosen by considering both temporal
molecular clustering and conventional age distributions. In `temporalCPSA`, we
implement a data-driven diagnostic approach: candidate classes are constrained
to ordered age intervals, while the AD-TMP matrix determines which intervals are
best supported.

These diagnostics are intended for cohort exploration and interpretation, not
as a required or decisive step in the CPSA pipeline.

The current beta exposes the underlying trajectory and clustering building
blocks while the full Figure 2 age-class reproduction workflow is being
assembled:

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

To review the ordered segmentation landscape, draw the green diagnostic heatmap.
The optimizer-selected solution can be highlighted, and manuscript-style or
investigator-selected class boundaries can be drawn on the same display. This
keeps the diagnostic honest: the green heatmap shows the molecular segmentation
landscape, while the overlaid cutpoints show the proposed age-class
interpretation.

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
predefined cutpoints after molecular age classes have been established, but they
are not a substitute for the TMP-based age-class exploration workflow.

## Manuscript Reproduction

The current beta also supports direct manuscript reproduction from public source
files:

Paper figure-generation scripts are maintained in the companion manuscript
repository:
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

The package is being organized around these modules:

- public data loading and sample harmonization from manuscript source tables;
- temporal molecular profile and trajectory generation;
- tumor versus normal/reference trajectory comparison, including Figure 2-style analyses;
- trajectory divergence summaries and visualization;
- clustering of age-dependent trajectory patterns;
- downstream association and survival modeling, including Cross-Population
  Survival Analysis (CPSA).

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

Normal/reference DLPFC developmental data used for trajectory analyses are documented package extdata derived from the DEveLopmental Trajectory Atlas (DELTA), PMID: 30518843, available at <http://amp.pharm.mssm.edu/DELTA>.

## Reproducibility Principle

Paper-facing workflows should read directly from public source files in `data/` whenever possible and avoid hidden `.RData` objects or manually generated intermediate TSV files. Serialized package data are reserved for documented external reference data that are part of the analysis provenance.

Detailed manuscript-specific trajectory settings are documented separately in
[`docs/manuscript-trajectory-reproduction-notes.md`](../docs/manuscript-trajectory-reproduction-notes.md).

The CPSA model engine computes raw p-values and, by default, Benjamini-Yekutieli
FDR columns with `stats::p.adjust(..., method = "BY")`. In the manuscript
reference-cohort analyses, protein and RNA significance calls used a separate
local-FDR procedure on transformed p-values; this is noted for reproducibility
context and is not built into the package API. The local-FDR specifications are
kept in the relevant manuscript figure/reproduction code.
