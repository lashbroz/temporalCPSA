# Package index

## Package overview

- [`temporalCPSA`](https://lashbroz.github.io/temporalCPSA/reference/temporalCPSA-package.md)
  [`temporalCPSA-package`](https://lashbroz.github.io/temporalCPSA/reference/temporalCPSA-package.md)
  : temporalCPSA: Cross-Population Survival Analysis via Temporal
  Molecular Profiles
- [`ageTMP_status()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_status.md)
  : Report temporalCPSA package status

## CPSA model specification and fitting

- [`ageTMP_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_cpsa_spec.md)
  : Create a CPSA model specification
- [`ageTMP_prepare_cdisc_cpsa_survival()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_prepare_cdisc_cpsa_survival.md)
  : Prepare reference-cohort clinical covariates for CPSA survival
  modeling
- [`ageTMP_prepare_reference_cpsa_survival()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_prepare_reference_cpsa_survival.md)
  : Prepare a reference-cohort CPSA survival model frame
- [`ageTMP_fit_reference_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_cpsa.md)
  : Fit reference-cohort CPSA Cox models
- [`ageTMP_fit_reference_protein_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_reference_protein_cpsa.md)
  : Fit reference-cohort protein CPSA models
- [`ageTMP_fit_cdisc_protein_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_cdisc_protein_cpsa.md)
  : Fit reference-cohort protein CPSA Cox models

## Temporal molecular profiles

- [`ageTMP_fit_tumor_trajectory()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_tumor_trajectory.md)
  : Fit tumor age-dependent molecular trajectories
- [`ageTMP_predict_tumor_trajectory_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_predict_tumor_trajectory_matrix.md)
  : Predict tumor age trajectories as a feature matrix
- [`ageTMP_compare_normal_tumor_trajectory()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_compare_normal_tumor_trajectory.md)
  : Compare normal and tumor age trajectories
- [`ageTMP_rank_trajectory_sd()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_rank_trajectory_sd.md)
  : Rank trajectory features by standard deviation
- [`ageTMP_filter_trajectory_sd()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_filter_trajectory_sd.md)
  : Filter trajectory features by standard deviation
- [`ageTMP_select_dynamic_tmp_features()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_select_dynamic_tmp_features.md)
  : Select dynamic TMP features

## Age-class exploration

- [`ageTMP_segment_age_classes()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_segment_age_classes.md)
  : Select contiguous age classes from an AD-TMP matrix
- [`ageTMP_segment_diagnostic_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_segment_diagnostic_matrix.md)
  : Build an ordered segment diagnostic matrix
- [`ageTMP_plot_segment_diagnostic()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_plot_segment_diagnostic.md)
  : Plot an ordered age-segmentation diagnostic heatmap
- [`ageTMP_derive_age_class()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_derive_age_class.md)
  : Derive age-class labels from age
- [`ageTMP_derive_age_class_range()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_derive_age_class_range.md)
  : Derive age-class interval labels from age
- [`ageTMP_add_age_class()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_age_class.md)
  : Add age-class columns to a clinical data frame
- [`ageTMP_combine_tmp_matrices()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_combine_tmp_matrices.md)
  : Combine TMP matrices for multi-omic sample clustering
- [`ageTMP_cluster_tmp_samples()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_cluster_tmp_samples.md)
  : Cluster samples from temporal molecular profiles
- [`ageTMP_summarize_tmp_age_clusters()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_summarize_tmp_age_clusters.md)
  : Summarize age cutpoints from TMP sample clusters

## Data loading and harmonization

- [`ageTMP_load_cdisc_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_clinical.md)
  : Load manuscript reference-cohort clinical data from STable1
- [`ageTMP_load_cdisc_glyco_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_glyco_matrix.md)
  : Load public cDisc glycopeptide data
- [`ageTMP_load_cdisc_mutation()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_mutation.md)
  : Load manuscript reference-cohort mutation calls for CPSA modeling
- [`ageTMP_load_cdisc_protein_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_protein_matrix.md)
  : Load manuscript reference-cohort protein data as a gene-by-sample
  matrix
- [`ageTMP_load_cdisc_rna_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_cdisc_rna_matrix.md)
  : Load cDisc RNA matrix
- [`ageTMP_load_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_clinical.md)
  : Load public manuscript clinical data
- [`ageTMP_load_discovery_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_discovery_clinical.md)
  : Load discovery-cohort clinical data
- [`ageTMP_load_feature_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_feature_matrix.md)
  : Load a public feature matrix
- [`ageTMP_load_figure1f_ad_tmp_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_figure1f_ad_tmp_matrix.md)
  : Load Figure 1F archived AD-TMP heatmap matrix
- [`ageTMP_load_figure1f_survival_annotation()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_figure1f_survival_annotation.md)
  : Load Figure 1F survival-days annotation data
- [`ageTMP_load_glyco_membrane_features()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_glyco_membrane_features.md)
  : Load membrane glycopeptide annotation
- [`ageTMP_load_molecular()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_molecular.md)
  : Load a public molecular data table
- [`ageTMP_load_normal_reference()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_normal_reference.md)
  : Load normal/reference molecular data used for AD-TMP trajectory
  analysis
- [`ageTMP_load_protein_trajectory_features()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_protein_trajectory_features.md)
  : Load manuscript protein trajectory feature universe
- [`ageTMP_load_reference_clinical()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_reference_clinical.md)
  : Load reference-cohort clinical data
- [`ageTMP_load_reference_mutation()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_reference_mutation.md)
  : Load reference-cohort mutation covariates
- [`ageTMP_load_reference_protein_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_reference_protein_matrix.md)
  : Load reference-cohort protein matrix
- [`ageTMP_load_supplement()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_supplement.md)
  : Read a public supplementary workbook sheet
- [`ageTMP_load_trajectory_qc()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_load_trajectory_qc.md)
  : Load manuscript trajectory QC omit lists
- [`ageTMP_data_sources()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_data_sources.md)
  : List public ageTMP manuscript data sources
- [`ageTMP_normalize_sample_ids()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_normalize_sample_ids.md)
  : Normalize HOPE AYA sample identifiers
- [`ageTMP_split_annotation_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_split_annotation_matrix.md)
  : Split a molecular table into annotation and numeric matrix
  components
- [`ageTMP_collapse_matrix_by_feature()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_collapse_matrix_by_feature.md)
  : Collapse matrix rows by feature identifier

## Manuscript reproduction helpers

- [`ageTMP_build_age_class_supplement()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_age_class_supplement.md)
  : Build the supplementary age-class table
- [`ageTMP_build_cdisc_glyco_adjusted_matrix()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_cdisc_glyco_adjusted_matrix.md)
  : Build source-derived protein-adjusted glycopeptide matrix
- [`ageTMP_build_sa_discovery_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_sa_discovery_cpsa.md)
  : Build manuscript STable4 discovery-cohort CPSA columns
- [`ageTMP_build_sa_glyco_disc_from_results()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_sa_glyco_disc_from_results.md)
  : Assemble manuscript STable4 SA-Glyco-Disc columns from glyco CPSA
  results
- [`ageTMP_build_sa_protein_cdisc()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_sa_protein_cdisc.md)
  : Build manuscript reference-cohort protein CPSA columns for STable4
  comparison
- [`ageTMP_build_sa_protein_reference()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_build_sa_protein_reference.md)
  : Build manuscript STable4 protein reference-cohort CPSA columns
- [`ageTMP_add_glyco_locfdr_by_manuscript_settings()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_glyco_locfdr_by_manuscript_settings.md)
  : Add glyco local-FDR columns with manuscript sex-specific settings
- [`ageTMP_add_glyco_locfdr_columns()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_glyco_locfdr_columns.md)
  : Add manuscript glyco local-FDR columns to a validation result
- [`ageTMP_add_glyco_ped_columns()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_add_glyco_ped_columns.md)
  : Add manuscript PED columns to a glyco ADO CPSA result
- [`ageTMP_prepare_glyco_ado_discovery_survival()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_prepare_glyco_ado_discovery_survival.md)
  : Prepare discovery glyco ADO-span survival data
- [`ageTMP_prepare_glyco_ado_reference_survival()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_prepare_glyco_ado_reference_survival.md)
  : Prepare reference glyco ADO-span survival data
- [`ageTMP_fit_glyco_ado_cpsa()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_glyco_ado_cpsa.md)
  : Fit manuscript glyco ADO-span CPSA models
- [`ageTMP_fit_glyco_ado_cpsa_with_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_fit_glyco_ado_cpsa_with_spec.md)
  : Fit glyco ADO CPSA models with an explicit specification
- [`ageTMP_glyco_ado_cpsa_spec()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_glyco_ado_cpsa_spec.md)
  : Manuscript glyco ADO Cox model specification
- [`ageTMP_glyco_locfdr_settings()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_glyco_locfdr_settings.md)
  : Return manuscript glyco local-FDR settings
- [`ageTMP_extract_glyco_trajectory_terms()`](https://lashbroz.github.io/temporalCPSA/reference/ageTMP_extract_glyco_trajectory_terms.md)
  : Extract sex-specific glycopeptide trajectory score matrices
