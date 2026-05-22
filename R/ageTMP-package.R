#' temporalCPSA: Cross-Population Survival Analysis via Temporal Molecular Profiles
#'
#' `temporalCPSA` provides tools for both manuscript reproduction and reusable
#' cross-population survival analysis via temporal molecular profiles.
#' The package scope includes TMP generation, normal/reference versus tumor
#' trajectory comparison, trajectory visualization and divergence analysis,
#' multi-omic age-class clustering, and downstream association or survival
#' modeling.
#'
#' @section Getting started:
#' Start with [ageTMP_cpsa_spec()] to define a survival model specification and
#' [ageTMP_fit_reference_cpsa()] to fit feature-wise Cox models in a reference
#' cohort. Users may supply their own age classes, use biologically motivated
#' clinical age classes, explore age classes from AD-TMP diagnostics, or omit
#' age-class structure entirely.
#'
#' @section Core functions:
#' \itemize{
#'   \item [ageTMP_cpsa_spec()] defines a CPSA model.
#'   \item [ageTMP_fit_reference_cpsa()] fits feature-wise reference-cohort
#'     survival models.
#'   \item [ageTMP_predict_tumor_trajectory_matrix()] derives temporal molecular
#'     profile matrices for downstream modeling.
#'   \item [ageTMP_rank_trajectory_sd()] and [ageTMP_filter_trajectory_sd()]
#'     screen low-dynamic-range trajectory features.
#'   \item [ageTMP_segment_age_classes()] and
#'     [ageTMP_plot_segment_diagnostic()] provide optional exploratory
#'     age-class diagnostics.
#' }
#'
#' @section Installed quick start:
#' A short installed Markdown guide is available at:
#' `system.file("doc", "temporalCPSA-quick-start.md", package = "temporalCPSA")`.
#' The GitHub README is available at
#' <https://github.com/lashbroz/temporalCPSA>.
#'
#' @examples
#' help(package = "temporalCPSA")
#' ?temporalCPSA
#' ?ageTMP_fit_reference_cpsa
"_PACKAGE"
