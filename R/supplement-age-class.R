#' Build the supplementary age-class table
#'
#' Build the `cDisc_AgeClass`-style supplementary table from clinical ages and,
#' when available, the TMP matrix used for silhouette scoring. Age-class labels
#' and ranges are derived from numeric ages. `sil_width` is computed only when a
#' feature-by-sample TMP matrix is supplied.
#'
#' @param clinical Clinical data frame.
#' @param tmp_matrix Optional feature-by-sample TMP matrix used to compute
#'   sample silhouette widths.
#' @param sample_ids Optional sample IDs to include. If `NULL`, all clinical
#'   samples with non-missing age are used, or the intersection with
#'   `colnames(tmp_matrix)` when a TMP matrix is supplied.
#' @param id_col Clinical sample ID column.
#' @param age_col Clinical age column.
#' @inheritParams ageTMP_derive_age_class
#'
#' @return A data frame with `id`, `age.class`, `age.class.range`, and
#'   `sil_width`.
#' @export
ageTMP_build_age_class_supplement <- function(
  clinical,
  tmp_matrix = NULL,
  sample_ids = NULL,
  id_col = "id",
  age_col = "age",
  breaks = c(0, 15, 26, 40, 62, 80),
  labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
) {
  clinical <- data.frame(clinical, stringsAsFactors = FALSE)
  required <- c(id_col, age_col)
  missing <- setdiff(required, names(clinical))
  if (length(missing) > 0) {
    stop("`clinical` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  ids <- ageTMP_normalize_sample_ids(clinical[[id_col]])
  age <- suppressWarnings(as.numeric(clinical[[age_col]]))
  names(age) <- ids
  age <- age[!is.na(age)]

  if (is.null(sample_ids)) {
    sample_ids <- names(age)
    if (!is.null(tmp_matrix)) {
      sample_ids <- intersect(ageTMP_normalize_sample_ids(colnames(tmp_matrix)), sample_ids)
    }
  } else {
    sample_ids <- intersect(ageTMP_normalize_sample_ids(sample_ids), names(age))
  }
  sample_ids <- sample_ids[order(age[sample_ids])]
  if (length(sample_ids) == 0) {
    stop("No samples remain after matching IDs and non-missing ages.", call. = FALSE)
  }

  out <- data.frame(
    id = sample_ids,
    age.class = as.character(ageTMP_derive_age_class(
      age[sample_ids],
      breaks = breaks,
      labels = labels,
      right = right,
      include_lowest = include_lowest,
      ordered = ordered
    )),
    age.class.range = as.character(ageTMP_derive_age_class_range(
      age[sample_ids],
      breaks = breaks,
      right = right,
      include_lowest = include_lowest,
      ordered = ordered
    )),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out$sil_width <- NA_real_
  if (!is.null(tmp_matrix)) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Package `cluster` is required to compute silhouette widths.", call. = FALSE)
    }
    tmp_matrix <- as.matrix(tmp_matrix)
    colnames(tmp_matrix) <- ageTMP_normalize_sample_ids(colnames(tmp_matrix))
    tmp_matrix <- tmp_matrix[, sample_ids, drop = FALSE]
    tmp_matrix <- tmp_matrix[rowSums(is.na(tmp_matrix)) == 0, colSums(is.na(tmp_matrix)) == 0, drop = FALSE]
    sample_ids2 <- colnames(tmp_matrix)
    classes <- out$age.class[match(sample_ids2, out$id)]
    if (length(unique(classes)) < 2) {
      stop("At least two age classes are required to compute silhouette widths.", call. = FALSE)
    }
    sil <- cluster::silhouette(as.numeric(factor(classes, levels = labels)), stats::dist(t(tmp_matrix)))
    out$sil_width[match(sample_ids2, out$id)] <- sil[, "sil_width"]
  }

  out
}
