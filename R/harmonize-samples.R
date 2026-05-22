#' Normalize HOPE AYA sample identifiers
#'
#' This mirrors the repeated `get_id()` helper used in the paper scripts.
#'
#' @param x Character vector of sample identifiers.
#'
#' @return A character vector with common prefixes removed and periods converted to dashes.
#' @export
ageTMP_normalize_sample_ids <- function(x) {
  x <- as.character(x)
  gsub("\\.", "-", gsub("^((X|A\\.|G\\.|P\\.))", "", x))
}

#' Split a molecular table into annotation and numeric matrix components
#'
#' @param data Molecular data frame with feature annotation columns followed by sample columns.
#' @param annotation_cols Integer or character vector identifying annotation columns.
#' @param row_id Optional annotation column to use as matrix row names.
#' @param normalize_colnames Whether to normalize sample IDs in matrix column names.
#'
#' @return A list with `annotation` and `matrix`.
#' @export
ageTMP_split_annotation_matrix <- function(
  data,
  annotation_cols,
  row_id = NULL,
  normalize_colnames = TRUE
) {
  if (missing(data) || missing(annotation_cols)) {
    stop("Both `data` and `annotation_cols` are required.", call. = FALSE)
  }

  annotation <- data[, annotation_cols, drop = FALSE]
  matrix_data <- data[, setdiff(seq_along(data), match(names(annotation), names(data))), drop = FALSE]
  matrix_data <- as.data.frame(lapply(matrix_data, function(x) suppressWarnings(as.numeric(x))))
  matrix_data <- as.matrix(matrix_data)

  if (!is.null(row_id)) {
    if (!row_id %in% names(annotation)) {
      stop("`row_id` not found in annotation columns: ", row_id, call. = FALSE)
    }
    rownames(matrix_data) <- annotation[[row_id]]
  }

  if (normalize_colnames) {
    colnames(matrix_data) <- ageTMP_normalize_sample_ids(colnames(matrix_data))
  }

  list(annotation = annotation, matrix = matrix_data)
}
