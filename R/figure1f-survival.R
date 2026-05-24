#' Load Figure 1F survival-days annotation data
#'
#' @description
#' Reads an external archived Figure 1F survival-days annotation TSV file. This
#' helper exists so reproduction scripts can load the Figure 1F descriptive
#' annotation explicitly, without treating it as a general CPSA survival
#' modeling input or bundling the manuscript-specific object as package data.
#'
#' When `quiet = FALSE`, the loader prints a short note explaining that this
#' object reflects the manuscript figure-preparation cohort rather than a final
#' source-derived survival modeling table.
#'
#' @param path Path to `figure1f_survival_annotation_data.tsv`.
#' @param quiet Logical; if `TRUE`, suppress the contextual message printed
#'   when the packaged annotation is loaded.
#' @return A data frame with the Figure 1F survival-days annotation data.
#' @export
ageTMP_load_figure1f_survival_annotation <- function(path, quiet = FALSE) {
  if (missing(path) || !nzchar(path) || !file.exists(path)) {
    stop("`path` must point to an existing Figure 1F survival-days annotation TSV file.", call. = FALSE)
  }
  if (!isTRUE(quiet)) {
    message(
      "Loading Figure 1F survival-days annotation data from `", path, "`. ",
      "This file reflects the reference/validation clinical cohort used ",
      "during manuscript figure preparation, including additional external ",
      "reference cases from CPTAC GBM studies and excluding study-associated ",
      "cases. It was used only as a descriptive survival-days annotation for ",
      "the Figure 1F visualization, not as a formal CPSA survival modeling input."
    )
  }
  utils::read.delim(path, check.names = FALSE)
}

#' Load Figure 1F archived AD-TMP heatmap matrix
#'
#' @description
#' Reads an external archived Figure 1F AD-TMP heatmap matrix from a TSV file.
#' This helper exists so manuscript reproduction scripts can load the archived
#' Figure 1F heatmap input with an explicit contextual note, while keeping the
#' matrix as a manuscript-specific reproduction input rather than bundled
#' package data.
#'
#' When `quiet = FALSE`, the loader prints a short note explaining that this
#' object reflects the manuscript figure-preparation matrix rather than a
#' source-derived matrix regenerated from final repository data tables.
#'
#' @param path Path to `figure1f_ad_tmp_legacy_clustme.tsv`, with feature names
#'   in the first column and sample IDs in the remaining columns.
#' @param quiet Logical; if `TRUE`, suppress the contextual message printed
#'   when the matrix is loaded.
#' @return A numeric matrix with AD-TMP features in rows and samples in columns.
#' @export
ageTMP_load_figure1f_ad_tmp_matrix <- function(path, quiet = FALSE) {
  if (missing(path) || !nzchar(path) || !file.exists(path)) {
    stop("`path` must point to an existing Figure 1F archived AD-TMP matrix TSV file.", call. = FALSE)
  }
  if (!isTRUE(quiet)) {
    message(
      "Loading Figure 1F archived AD-TMP heatmap matrix from `", path, "`. ",
      "This file reflects the manuscript figure-preparation input and is used ",
      "for exact visualization reproduction. For a final-data analogue, use ",
      "the source-derived Figure 1F script."
    )
  }
  mat <- utils::read.delim(path, check.names = FALSE)
  rownames(mat) <- mat[[1]]
  mat <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}
