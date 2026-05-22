#' Collapse matrix rows by feature identifier
#'
#' @param mat Numeric feature-by-sample matrix.
#' @param feature Character vector assigning each row of `mat` to a collapsed
#'   feature identifier.
#' @param method Summary method. Currently only `"mean"` is supported.
#'
#' @return A numeric matrix with one row per unique feature.
#' @export
ageTMP_collapse_matrix_by_feature <- function(mat, feature, method = "mean") {
  method <- match.arg(method, "mean")
  mat <- as.matrix(mat)
  feature <- as.character(feature)

  if (nrow(mat) != length(feature)) {
    stop("`feature` must have one value per row of `mat`.", call. = FALSE)
  }

  keep <- !is.na(feature) & nzchar(feature)
  mat <- mat[keep, , drop = FALSE]
  feature <- feature[keep]

  collapsed <- lapply(split(seq_along(feature), feature), function(i) {
    colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
  })

  out <- do.call(rbind, collapsed)
  rownames(out) <- names(collapsed)
  out[is.nan(out)] <- NA_real_
  out
}
