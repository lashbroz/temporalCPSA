select_adaptive_span <- function(x, y, min_span = 0.5, max_span = 3, span_step = 0.1) {
  if (!requireNamespace("locfit", quietly = TRUE)) {
    stop("Package `locfit` is required when `adaptive_span = TRUE`.", call. = FALSE)
  }

  ok <- !is.na(x) & !is.na(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  if (length(unique(x)) < 4 || length(y) < 6) {
    return(NA_real_)
  }

  spans <- seq(min_span, max_span, by = span_step)
  data <- data.frame(x = x, y = y)
  scores <- vapply(spans, function(alpha) {
    score <- try(locfit::gcv(y ~ x, data = data, alpha = alpha), silent = TRUE)
    if (inherits(score, "try-error")) {
      return(Inf)
    }
    round(as.numeric(score["gcv"]), 3)
  }, numeric(1))

  if (all(!is.finite(scores))) {
    return(NA_real_)
  }
  spans[which.min(scores)]
}

#' Fit tumor age-dependent molecular trajectories
#'
#' Fit sex-stratified age-dependent tumor molecular trajectories from a public
#' feature-by-sample matrix and sample metadata. This is the tumor-only AD-TMP
#' fitting step used by manuscript heatmap-style trajectory panels such as
#' Figure 2A.
#'
#' @details
#' The function intentionally mirrors key details from the manuscript tumor
#' trajectory workflow. Tumor sample IDs are harmonized with
#' [ageTMP_normalize_sample_ids()], rows are centered/scaled using the requested
#' `center_age_range`, optional fitting can be restricted with `fit_age_range`,
#' and sex-stratified [stats::loess()] models are fit with base loess defaults.
#'
#' For exact manuscript reproduction, the original protein workflow used
#' feature/sex-specific adaptive loess spans selected by generalized
#' cross-validation over a span grid. Set `adaptive_span = TRUE` to recompute
#' those spans from the public tumor matrix; this requires the suggested
#' `locfit` package and can be slow for thousands of proteins. A single numeric
#' `span`, or a data frame with `feature`, `sex`, `tissue`, and `span` columns,
#' can also be supplied when spans are known or when a faster reconstruction is
#' desired.
#'
#' @param tumor_mat Tumor feature-by-sample matrix.
#' @param tumor_metadata Tumor sample metadata.
#' @param features Features to model. Defaults to all matrix rows.
#' @param tumor_sample_col Tumor metadata sample ID column.
#' @param tumor_age_col Tumor metadata age column.
#' @param tumor_sex_col Tumor metadata sex column.
#' @param center_age_range Age range used to center and scale each feature.
#' @param fit_age_range Optional age range used for loess fitting.
#' @param pre_scale Whether to row-center and row-scale the tumor matrix before
#'   the age-range centering step. This mirrors manuscript trajectory scripts
#'   that first standardized the complete feature matrix, then standardized
#'   again within the modeled age range.
#' @param span Loess span. May be a single number or a span data frame accepted
#'   by the normal/tumor trajectory functions.
#' @param adaptive_span Recompute feature/sex-specific spans by GCV.
#' @param min_span Minimum span for adaptive selection.
#' @param max_span Maximum span for adaptive selection.
#' @param span_step Span grid step for adaptive selection.
#' @param prediction_ages Optional ages where trajectories are predicted. If
#'   `NULL`, predictions are made at all harmonized tumor sample ages.
#' @param prediction_sample_ids Optional IDs corresponding to `prediction_ages`.
#' @param ci_level Confidence level for fitted trajectories.
#' @param n_cores Number of parallel worker processes for feature/sex fits on
#'   Unix-like systems. Use `1` for serial execution.
#' @param progress Whether to print simple progress messages.
#'
#' @return A data frame with one row per feature, sex, and prediction age.
#' @export
ageTMP_fit_tumor_trajectory <- function(
  tumor_mat,
  tumor_metadata,
  features = rownames(tumor_mat),
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  center_age_range = c(0, 26),
  fit_age_range = NULL,
  pre_scale = FALSE,
  span = 1.5,
  adaptive_span = FALSE,
  min_span = 0.5,
  max_span = 3,
  span_step = 0.1,
  prediction_ages = NULL,
  prediction_sample_ids = NULL,
  ci_level = 0.95,
  n_cores = 1,
  progress = FALSE
) {
  tumor_mat <- as.matrix(tumor_mat)
  features <- intersect(features, rownames(tumor_mat))
  if (length(features) == 0) {
    stop("No requested features are present in `tumor_mat`.", call. = FALSE)
  }

  required_cols <- c(tumor_sample_col, tumor_age_col, tumor_sex_col)
  missing_cols <- setdiff(required_cols, names(tumor_metadata))
  if (length(missing_cols) > 0) {
    stop(
      "Tumor metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  tumor_metadata <- data.frame(tumor_metadata, stringsAsFactors = FALSE)
  tumor_metadata[[tumor_sample_col]] <- ageTMP_normalize_sample_ids(tumor_metadata[[tumor_sample_col]])
  tumor_metadata[[tumor_sex_col]] <- standardize_sex(tumor_metadata[[tumor_sex_col]])
  tumor_metadata[[tumor_age_col]] <- as.numeric(tumor_metadata[[tumor_age_col]])
  colnames(tumor_mat) <- ageTMP_normalize_sample_ids(colnames(tumor_mat))

  tumor_ids <- intersect(colnames(tumor_mat), tumor_metadata[[tumor_sample_col]])
  if (length(tumor_ids) == 0) {
    stop("No tumor matrix columns match tumor metadata sample IDs.", call. = FALSE)
  }

  tumor_metadata <- tumor_metadata[match(tumor_ids, tumor_metadata[[tumor_sample_col]]), , drop = FALSE]
  tumor_mat <- tumor_mat[features, tumor_ids, drop = FALSE]

  fit_keep <- !is.na(tumor_metadata[[tumor_age_col]]) & tumor_metadata[[tumor_sex_col]] %in% c("Male", "Female")
  if (!is.null(fit_age_range)) {
    fit_keep <- fit_keep &
      tumor_metadata[[tumor_age_col]] >= fit_age_range[1] &
      tumor_metadata[[tumor_age_col]] <= fit_age_range[2]
  }
  if (!any(fit_keep)) {
    stop("No tumor samples remain after age/sex filtering.", call. = FALSE)
  }

  fit_metadata <- tumor_metadata[fit_keep, , drop = FALSE]
  fit_mat <- tumor_mat[, fit_keep, drop = FALSE]
  if (isTRUE(pre_scale)) {
    fit_mat <- t(scale(t(fit_mat)))
  }

  center_keep <- fit_metadata[[tumor_age_col]] >= center_age_range[1] &
    fit_metadata[[tumor_age_col]] <= center_age_range[2]
  center_ids <- fit_metadata[[tumor_sample_col]][center_keep]
  tumor_scaled <- scale_rows(fit_mat, center_ids = center_ids)

  if (is.null(prediction_ages)) {
    pred_age <- tumor_metadata[[tumor_age_col]]
    pred_sample_id <- tumor_metadata[[tumor_sample_col]]
  } else {
    pred_age <- as.numeric(prediction_ages)
    if (is.null(prediction_sample_ids)) {
      pred_sample_id <- rep(NA_character_, length(pred_age))
    } else {
      if (length(prediction_sample_ids) != length(pred_age)) {
        stop("`prediction_sample_ids` must have the same length as `prediction_ages`.", call. = FALSE)
      }
      pred_sample_id <- ageTMP_normalize_sample_ids(prediction_sample_ids)
    }
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  pairs <- expand.grid(feature = features, sex = c("Male", "Female"), stringsAsFactors = FALSE)
  fit_pair <- function(i) {
      feature <- pairs$feature[[i]]
      sex <- pairs$sex[[i]]
      if (isTRUE(progress) && (i == 1L || i %% 500L == 0L || i == nrow(pairs))) {
        message("ageTMP trajectory fit ", i, "/", nrow(pairs))
      }
      sex_idx <- which(fit_metadata[[tumor_sex_col]] == sex)
      if (length(sex_idx) == 0) {
        return(NULL)
      }
      tumor_age <- fit_metadata[[tumor_age_col]][sex_idx]
      tumor_values <- as.numeric(tumor_scaled[feature, sex_idx])

      tumor_span <- if (isTRUE(adaptive_span)) {
        select_adaptive_span(tumor_age, tumor_values, min_span, max_span, span_step)
      } else {
        resolve_trajectory_span(span, feature, sex, "Tumor")
      }

      pred <- if (is.na(tumor_span)) {
        list(fit = rep(NA_real_, length(pred_age)), se = rep(NA_real_, length(pred_age)))
      } else {
        fit_loess_predict(tumor_age, tumor_values, pred_age, span = tumor_span)
      }

      data.frame(
        feature = feature,
        sample_id = pred_sample_id,
        age = pred_age,
        sex = sex,
        tissue = "Tumor",
        fit = pred$fit,
        se = pred$se,
        ci_lower = pred$fit - z * pred$se,
        ci_upper = pred$fit + z * pred$se,
        span = tumor_span,
        stringsAsFactors = FALSE
      )
  }

  n_cores <- max(1L, as.integer(n_cores)[1])
  out <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(pairs)), fit_pair, mc.cores = n_cores)
  } else {
    lapply(seq_len(nrow(pairs)), fit_pair)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  do.call(rbind, out)
}

#' Predict tumor age trajectories as a feature matrix
#'
#' This helper wraps [ageTMP_fit_tumor_trajectory()] for CPSA workflows where
#' the age-dependent tumor trajectory is evaluated at the ages of an external
#' reference cohort. The result is a numeric feature-by-sample matrix that can be
#' passed directly to [ageTMP_fit_reference_cpsa()].
#'
#' @inheritParams ageTMP_fit_tumor_trajectory
#' @param prediction_metadata Data frame containing sample IDs, ages, and sex
#'   labels for the cohort where trajectory values should be predicted.
#' @param prediction_sample_col Sample ID column in `prediction_metadata`.
#' @param prediction_age_col Age column in `prediction_metadata`.
#' @param prediction_sex_col Sex column in `prediction_metadata`.
#' @param prediction_scope Whether each sex-specific trajectory should be
#'   predicted only for samples with the matching sex label (`"matching_sex"`)
#'   or for all samples in `prediction_metadata` (`"all_samples"`). The latter
#'   is useful when constructing sex-stratified trajectory rows over a common
#'   sample grid.
#' @param return_trajectory Whether to include the long-form trajectory table in
#'   the result. Set to `FALSE` for large CPSA matrices where only the
#'   feature-by-sample prediction matrix is needed.
#'
#' @return A named list with `matrix` and the long-form `trajectory` table.
#' @export
ageTMP_predict_tumor_trajectory_matrix <- function(
  tumor_mat,
  tumor_metadata,
  prediction_metadata,
  features = rownames(tumor_mat),
  tumor_sample_col = "id",
  tumor_age_col = "age",
  tumor_sex_col = "sex",
  prediction_sample_col = "id",
  prediction_age_col = "age",
  prediction_sex_col = "Gender",
  center_age_range = c(0, 62),
  fit_age_range = c(0, 80),
  pre_scale = FALSE,
  span = 1.5,
  adaptive_span = FALSE,
  min_span = 0.5,
  max_span = 2,
  span_step = 0.1,
  ci_level = 0.95,
  n_cores = 1,
  progress = FALSE,
  prediction_scope = c("matching_sex", "all_samples"),
  return_trajectory = TRUE
) {
  prediction_scope <- match.arg(prediction_scope)
  prediction_metadata <- data.frame(prediction_metadata, stringsAsFactors = FALSE)
  required <- c(prediction_sample_col, prediction_age_col, prediction_sex_col)
  missing <- setdiff(required, names(prediction_metadata))
  if (length(missing) > 0) {
    stop(
      "Prediction metadata is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  pred_ids <- ageTMP_normalize_sample_ids(prediction_metadata[[prediction_sample_col]])
  pred_age <- suppressWarnings(as.numeric(prediction_metadata[[prediction_age_col]]))
  pred_sex <- standardize_sex(prediction_metadata[[prediction_sex_col]])

  if (!isTRUE(return_trajectory)) {
    tumor_mat <- as.matrix(tumor_mat)
    features <- intersect(features, rownames(tumor_mat))
    if (length(features) == 0) {
      stop("No requested features are present in `tumor_mat`.", call. = FALSE)
    }

    tumor_metadata <- data.frame(tumor_metadata, stringsAsFactors = FALSE)
    tumor_metadata[[tumor_sample_col]] <- ageTMP_normalize_sample_ids(tumor_metadata[[tumor_sample_col]])
    tumor_metadata[[tumor_sex_col]] <- standardize_sex(tumor_metadata[[tumor_sex_col]])
    tumor_metadata[[tumor_age_col]] <- as.numeric(tumor_metadata[[tumor_age_col]])
    colnames(tumor_mat) <- ageTMP_normalize_sample_ids(colnames(tumor_mat))

    tumor_ids <- intersect(colnames(tumor_mat), tumor_metadata[[tumor_sample_col]])
    if (length(tumor_ids) == 0) {
      stop("No tumor matrix columns match tumor metadata sample IDs.", call. = FALSE)
    }

    tumor_metadata <- tumor_metadata[match(tumor_ids, tumor_metadata[[tumor_sample_col]]), , drop = FALSE]
    tumor_mat <- tumor_mat[features, tumor_ids, drop = FALSE]

    fit_keep <- !is.na(tumor_metadata[[tumor_age_col]]) & tumor_metadata[[tumor_sex_col]] %in% c("Male", "Female")
    if (!is.null(fit_age_range)) {
      fit_keep <- fit_keep &
        tumor_metadata[[tumor_age_col]] >= fit_age_range[1] &
        tumor_metadata[[tumor_age_col]] <= fit_age_range[2]
    }
    if (!any(fit_keep)) {
      stop("No tumor samples remain after age/sex filtering.", call. = FALSE)
    }

    fit_metadata <- tumor_metadata[fit_keep, , drop = FALSE]
    fit_mat <- tumor_mat[, fit_keep, drop = FALSE]
    if (isTRUE(pre_scale)) {
      fit_mat <- t(scale(t(fit_mat)))
    }

    center_keep <- fit_metadata[[tumor_age_col]] >= center_age_range[1] &
      fit_metadata[[tumor_age_col]] <= center_age_range[2]
    center_ids <- fit_metadata[[tumor_sample_col]][center_keep]
    tumor_scaled <- scale_rows(fit_mat, center_ids = center_ids)

    keep <- !is.na(pred_ids) & !is.na(pred_age) & pred_sex %in% c("Male", "Female")
    pred_ids_keep <- pred_ids[keep]
    pred_age_keep <- pred_age[keep]
    pred_sex_keep <- pred_sex[keep]

    pairs <- expand.grid(feature = features, sex = c("Male", "Female"), stringsAsFactors = FALSE)
    fit_pair <- function(i) {
      feature <- pairs$feature[[i]]
      sex <- pairs$sex[[i]]
      if (isTRUE(progress) && (i == 1L || i %% 500L == 0L || i == nrow(pairs))) {
        message("ageTMP trajectory fit ", i, "/", nrow(pairs))
      }
      sex_idx <- which(fit_metadata[[tumor_sex_col]] == sex)
      pred_idx <- if (identical(prediction_scope, "all_samples")) {
        seq_along(pred_ids_keep)
      } else {
        which(pred_sex_keep == sex)
      }
      if (length(sex_idx) == 0 || length(pred_idx) == 0) {
        return(list(feature = feature, sex = sex, sample_id = character(), fit = numeric(), span = NA_real_))
      }
      tumor_age <- fit_metadata[[tumor_age_col]][sex_idx]
      tumor_values <- as.numeric(tumor_scaled[feature, sex_idx])

      tumor_span <- if (isTRUE(adaptive_span)) {
        select_adaptive_span(tumor_age, tumor_values, min_span, max_span, span_step)
      } else {
        resolve_trajectory_span(span, feature, sex, "Tumor")
      }
      fit <- if (is.na(tumor_span)) {
        rep(NA_real_, length(pred_idx))
      } else {
        fit_loess_predict(tumor_age, tumor_values, pred_age_keep[pred_idx], span = tumor_span)$fit
      }
      list(
        feature = feature,
        sex = sex,
        sample_id = pred_ids_keep[pred_idx],
        fit = fit,
        span = tumor_span
      )
    }

    n_cores <- max(1L, as.integer(n_cores)[1])
    fits <- if (n_cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(seq_len(nrow(pairs)), fit_pair, mc.cores = n_cores)
    } else {
      lapply(seq_len(nrow(pairs)), fit_pair)
    }

    out_rows <- if (identical(prediction_scope, "all_samples")) {
      paste(pairs$sex, pairs$feature, sep = "::")
    } else {
      features
    }
    out <- matrix(NA_real_, nrow = length(unique(out_rows)), ncol = length(pred_ids_keep))
    rownames(out) <- unique(out_rows)
    colnames(out) <- pred_ids_keep
    spans <- data.frame(
      feature = vapply(fits, `[[`, character(1), "feature"),
      sex = vapply(fits, `[[`, character(1), "sex"),
      span = vapply(fits, `[[`, numeric(1), "span"),
      stringsAsFactors = FALSE
    )
    for (fit in fits) {
      if (length(fit$sample_id) > 0) {
        out_row <- if (identical(prediction_scope, "all_samples")) {
          paste(fit$sex, fit$feature, sep = "::")
        } else {
          fit$feature
        }
        out[out_row, fit$sample_id] <- fit$fit
      }
    }
    return(list(matrix = out, trajectory = NULL, spans = spans))
  }

  trajectory <- ageTMP_fit_tumor_trajectory(
    tumor_mat = tumor_mat,
    tumor_metadata = tumor_metadata,
    features = features,
    tumor_sample_col = tumor_sample_col,
    tumor_age_col = tumor_age_col,
    tumor_sex_col = tumor_sex_col,
    center_age_range = center_age_range,
    fit_age_range = fit_age_range,
    pre_scale = pre_scale,
    span = span,
    adaptive_span = adaptive_span,
    min_span = min_span,
    max_span = max_span,
    span_step = span_step,
    prediction_ages = pred_age,
    prediction_sample_ids = pred_ids,
    ci_level = ci_level,
    n_cores = n_cores,
    progress = progress
  )

  keep <- !is.na(pred_ids) & !is.na(pred_age) & pred_sex %in% c("Male", "Female")
  pred_ids <- pred_ids[keep]
  pred_sex <- pred_sex[keep]

  feature_levels <- unique(trajectory$feature)
  out <- matrix(NA_real_, nrow = length(feature_levels), ncol = length(pred_ids))
  rownames(out) <- feature_levels
  colnames(out) <- pred_ids

  for (sex in c("Male", "Female")) {
    sample_ids <- pred_ids[pred_sex == sex]
    if (length(sample_ids) == 0) {
      next
    }
    sub <- trajectory[trajectory$sex == sex & trajectory$sample_id %in% sample_ids, , drop = FALSE]
    split_fit <- split(sub$fit, list(sub$feature, sub$sample_id), drop = TRUE)
    for (nm in names(split_fit)) {
      key <- strsplit(nm, ".", fixed = TRUE)[[1]]
      if (length(key) == 2 && key[[1]] %in% rownames(out) && key[[2]] %in% colnames(out)) {
        out[key[[1]], key[[2]]] <- split_fit[[nm]][1]
      }
    }
  }

  list(matrix = out, trajectory = trajectory)
}
