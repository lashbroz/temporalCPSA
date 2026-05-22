#' Select contiguous age classes from an AD-TMP matrix
#'
#' Partition age-ordered samples into contiguous age segments using molecular
#' structure in an age-dependent temporal molecular profile (AD-TMP) matrix.
#' Unlike unconstrained clustering, this procedure only considers segmentations
#' that can be written as ordered age intervals. It therefore guarantees
#' age-contiguous classes by construction.
#'
#' The method is a reusable formalization inspired by the manuscript
#' age-class workflow, where AD-TMP clustering and an age-ordered green
#' diagnostic heatmap were used to judge molecularly supported age strata. The
#' original manuscript step was interpretive: age-contiguous structure was
#' assessed visually from the ordered clustering diagnostics rather than solved
#' as a constrained optimization problem. In the manuscript, the final age
#' classes balanced data-driven AD-TMP structure with developmental precedent,
#' with particular attention to retaining adolescence as a distinct
#' biology-informed interval. This function makes the contiguity principle
#' explicit for new cohorts by optimizing over possible age cutpoints directly.
#' It is intended as an exploratory diagnostic for age-class assessment and
#' interpretation, not as a decisive automated step in the CPSA pipeline.
#'
#' @param tmp_matrix Numeric AD-TMP matrix with features in rows and samples in
#'   columns.
#' @param sample_age Numeric sample ages. If named, names are matched to
#'   `colnames(tmp_matrix)` after sample-ID normalization. If unnamed, the
#'   vector must be in the same order as `colnames(tmp_matrix)`.
#' @param k Number of contiguous age segments to select.
#' @param min_n Minimum number of samples allowed in each segment.
#' @param labels Optional labels for the resulting age classes. If `NULL`,
#'   labels are `segment1`, `segment2`, ...
#' @param scale_rows Whether to center and scale AD-TMP rows before scoring
#'   segments. Set this to `FALSE` when `tmp_matrix` has already been row-scaled.
#'
#' @return An object of class `ageTMP_age_segmentation`, a list containing the
#'   ordered sample metadata, named class assignments, selected cutpoints,
#'   segment summary, dynamic-programming score, and the fitted cost matrix.
#' @export
ageTMP_segment_age_classes <- function(
  tmp_matrix,
  sample_age,
  k = 5,
  min_n = 10,
  labels = NULL,
  scale_rows = FALSE
) {
  tmp_matrix <- as.matrix(tmp_matrix)
  if (is.null(rownames(tmp_matrix)) || is.null(colnames(tmp_matrix))) {
    stop("`tmp_matrix` must have row and column names.", call. = FALSE)
  }
  if (isTRUE(scale_rows)) {
    tmp_matrix <- t(scale(t(tmp_matrix)))
  }
  tmp_matrix <- tmp_matrix[rowSums(is.na(tmp_matrix)) == 0, , drop = FALSE]
  if (nrow(tmp_matrix) < 1) {
    stop("`tmp_matrix` must contain at least one complete row.", call. = FALSE)
  }

  k <- as.integer(k)
  min_n <- as.integer(min_n)
  if (length(k) != 1 || is.na(k) || k < 1) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  if (length(min_n) != 1 || is.na(min_n) || min_n < 1) {
    stop("`min_n` must be a positive integer.", call. = FALSE)
  }

  sample_age <- suppressWarnings(as.numeric(sample_age))
  if (is.null(names(sample_age))) {
    if (length(sample_age) != ncol(tmp_matrix)) {
      stop("Unnamed `sample_age` must have length `ncol(tmp_matrix)`.", call. = FALSE)
    }
    names(sample_age) <- colnames(tmp_matrix)
  } else {
    names(sample_age) <- ageTMP_normalize_sample_ids(names(sample_age))
  }
  matrix_ids <- ageTMP_normalize_sample_ids(colnames(tmp_matrix))
  age <- sample_age[match(matrix_ids, names(sample_age))]
  keep <- !is.na(age)
  tmp_matrix <- tmp_matrix[, keep, drop = FALSE]
  matrix_ids <- matrix_ids[keep]
  age <- age[keep]
  if (ncol(tmp_matrix) < k * min_n) {
    stop("Not enough samples for `k` segments with `min_n` samples each.", call. = FALSE)
  }

  ord <- order(age, matrix_ids)
  tmp_matrix <- tmp_matrix[, ord, drop = FALSE]
  matrix_ids <- matrix_ids[ord]
  age <- age[ord]
  n <- length(age)

  costs <- ageTMP_interval_wss_costs(t(tmp_matrix))
  dp <- matrix(Inf, nrow = k, ncol = n)
  prev <- matrix(NA_integer_, nrow = k, ncol = n)

  for (j in seq.int(min_n, n)) {
    dp[1, j] <- costs[1, j]
  }
  if (k > 1) {
    for (seg in seq.int(2, k)) {
      j_min <- seg * min_n
      for (j in seq.int(j_min, n)) {
        i_min <- (seg - 1L) * min_n
        i_max <- j - min_n
        candidates <- dp[seg - 1L, i_min:i_max] + costs[(i_min:i_max) + 1L, j]
        best <- which.min(candidates)
        dp[seg, j] <- candidates[best]
        prev[seg, j] <- (i_min:i_max)[best]
      }
    }
  }

  if (!is.finite(dp[k, n])) {
    stop("No valid age segmentation found.", call. = FALSE)
  }
  ends <- integer(k)
  starts <- integer(k)
  ends[k] <- n
  if (k > 1) {
    for (seg in seq.int(k, 2L)) {
      ends[seg - 1L] <- prev[seg, ends[seg]]
    }
  }
  starts[1] <- 1L
  if (k > 1) {
    starts[2:k] <- ends[seq_len(k - 1L)] + 1L
  }

  if (is.null(labels)) {
    labels <- paste0("segment", seq_len(k))
  }
  if (length(labels) != k) {
    stop("`labels` must have length `k`.", call. = FALSE)
  }

  segment_id <- integer(n)
  for (seg in seq_len(k)) {
    segment_id[starts[seg]:ends[seg]] <- seg
  }
  classes <- factor(labels[segment_id], levels = labels, ordered = TRUE)
  names(classes) <- matrix_ids

  summary <- data.frame(
    segment = seq_len(k),
    label = labels,
    start_index = starts,
    end_index = ends,
    n = ends - starts + 1L,
    min_age = age[starts],
    max_age = age[ends],
    stringsAsFactors = FALSE
  )
  cutpoints <- if (k > 1) age[ends[seq_len(k - 1L)]] else numeric()
  ordered_samples <- data.frame(
    id = matrix_ids,
    age = age,
    segment = segment_id,
    label = classes,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      classes = classes,
      cutpoints = cutpoints,
      segment_summary = summary,
      ordered_samples = ordered_samples,
      score = dp[k, n],
      k = k,
      min_n = min_n,
      cost_matrix = costs,
      dp = dp
    ),
    class = "ageTMP_age_segmentation"
  )
}

#' Build an ordered segment diagnostic matrix
#'
#' Fit contiguous AD-TMP segmentations for a range of `k` values and return a
#' matrix whose rows are `k` values and whose columns are age-ordered samples.
#' This provides the data needed to draw a green age-ordered diagnostic heatmap
#' analogous to the manuscript age-class selection display, while preserving
#' contiguity for every fitted row.
#'
#' @inheritParams ageTMP_segment_age_classes
#' @param k_values Integer vector of segment counts to fit.
#'
#' @return A list with `matrix`, `ordered_samples`, and per-`k` `fits`.
#' @export
ageTMP_segment_diagnostic_matrix <- function(
  tmp_matrix,
  sample_age,
  k_values = 2:6,
  min_n = 10,
  scale_rows = FALSE
) {
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[!is.na(k_values) & k_values >= 1]
  if (length(k_values) == 0) {
    stop("`k_values` must contain at least one positive integer.", call. = FALSE)
  }
  fits <- lapply(k_values, function(k) {
    ageTMP_segment_age_classes(
      tmp_matrix = tmp_matrix,
      sample_age = sample_age,
      k = k,
      min_n = min_n,
      scale_rows = scale_rows
    )
  })
  names(fits) <- paste0("k", k_values)
  ordered_ids <- fits[[1]]$ordered_samples$id
  mat <- do.call(rbind, lapply(fits, function(fit) {
    out <- fit$ordered_samples$segment
    names(out) <- fit$ordered_samples$id
    out[ordered_ids]
  }))
  rownames(mat) <- names(fits)
  colnames(mat) <- ordered_ids
  list(matrix = mat, ordered_samples = fits[[1]]$ordered_samples, fits = fits)
}

#' Plot an ordered age-segmentation diagnostic heatmap
#'
#' Draw the green ordered-segmentation diagnostic used to evaluate candidate
#' age-class resolutions. Rows show different values of `k`; columns are samples
#' sorted by age; colors encode contiguous segment membership. The selected
#' solution is overlaid with vertical cutpoint guides so the proposed age classes
#' can be read directly from the diagnostic display.
#'
#' @param diagnostic Output of [ageTMP_segment_diagnostic_matrix()].
#' @param selected_k Which fitted `k` to overlay. Defaults to the largest fitted
#'   `k`.
#' @param class_labels Optional labels for the selected classes. If `NULL`, the
#'   labels stored in the selected fit are used.
#' @param palette Optional color palette. Defaults to a light-to-dark green
#'   palette.
#' @param main Plot title.
#' @param show_age_axis Whether to show an age axis below the heatmap.
#' @param suggested_cutpoints Optional numeric age cutpoints to overlay across
#'   the diagnostic heatmap. Use this to draw proposed or manuscript-selected
#'   age classes on top of the ordered segmentation landscape.
#' @param suggested_labels Optional labels for the age intervals defined by
#'   `suggested_cutpoints`.
#' @param ... Additional arguments passed to [graphics::image()].
#'
#' @return Invisibly returns `diagnostic`.
#' @export
ageTMP_plot_segment_diagnostic <- function(
  diagnostic,
  selected_k = NULL,
  class_labels = NULL,
  palette = NULL,
  main = "AD-TMP contiguous age-segmentation diagnostic",
  show_age_axis = TRUE,
  suggested_cutpoints = NULL,
  suggested_labels = NULL,
  ...
) {
  if (!is.list(diagnostic) || !all(c("matrix", "ordered_samples", "fits") %in% names(diagnostic))) {
    stop("`diagnostic` must be the output of `ageTMP_segment_diagnostic_matrix()`.", call. = FALSE)
  }
  mat <- as.matrix(diagnostic$matrix)
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("`diagnostic$matrix` is empty.", call. = FALSE)
  }
  if (is.null(palette)) {
    palette <- grDevices::colorRampPalette(c("#d9f0d3", "#74c476", "#006d2c"))(max(mat, na.rm = TRUE))
  }
  if (is.null(selected_k)) {
    selected_name <- utils::tail(names(diagnostic$fits), 1)
  } else {
    selected_name <- paste0("k", as.integer(selected_k))
  }
  if (!selected_name %in% names(diagnostic$fits)) {
    stop("Selected `k` was not fit in `diagnostic`.", call. = FALSE)
  }
  selected <- diagnostic$fits[[selected_name]]
  age <- diagnostic$ordered_samples$age

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  graphics::par(mar = c(if (isTRUE(show_age_axis)) 4 else 2, 5, 4, 1))
  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = palette,
    axes = FALSE,
    xlab = if (isTRUE(show_age_axis)) "Age" else "",
    ylab = "Candidate k",
    main = main,
    ...
  )
  graphics::axis(
    side = 2,
    at = seq_len(nrow(mat)),
    labels = rev(rownames(mat)),
    las = 1
  )
  if (isTRUE(show_age_axis)) {
    ticks <- pretty(age)
    tick_pos <- vapply(ticks, function(x) which.min(abs(age - x)), integer(1))
    keep <- !duplicated(tick_pos) & ticks >= min(age, na.rm = TRUE) & ticks <= max(age, na.rm = TRUE)
    graphics::axis(side = 1, at = tick_pos[keep], labels = ticks[keep])
  }

  if (!is.null(suggested_cutpoints)) {
    suggested_cutpoints <- sort(unique(suppressWarnings(as.numeric(suggested_cutpoints))))
    suggested_cutpoints <- suggested_cutpoints[
      !is.na(suggested_cutpoints) &
        suggested_cutpoints > min(age, na.rm = TRUE) &
        suggested_cutpoints < max(age, na.rm = TRUE)
    ]
    suggested_pos <- vapply(suggested_cutpoints, function(x) {
      sum(age <= x) + 0.5
    }, numeric(1))
    if (length(suggested_pos) > 0) {
      graphics::abline(v = suggested_pos, col = "#111111", lwd = 2, lty = 1)
      graphics::mtext(
        text = format(suggested_cutpoints, trim = TRUE),
        side = 3,
        at = suggested_pos,
        line = 0.15,
        cex = 0.75
      )
    }
    if (!is.null(suggested_labels)) {
      breaks <- c(min(age, na.rm = TRUE), suggested_cutpoints, max(age, na.rm = TRUE))
      if (length(suggested_labels) != length(breaks) - 1L) {
        stop("`suggested_labels` must have one label per suggested age interval.", call. = FALSE)
      }
      interval_mids <- vapply(seq_len(length(breaks) - 1L), function(i) {
        idx <- which(age >= breaks[i] & age <= breaks[i + 1L])
        if (length(idx) == 0) {
          return(NA_real_)
        }
        mean(range(idx))
      }, numeric(1))
      graphics::mtext(
        text = suggested_labels[!is.na(interval_mids)],
        side = 3,
        at = interval_mids[!is.na(interval_mids)],
        line = 1.05,
        cex = 0.8,
        font = 2
      )
    }
  }

  selected_y <- nrow(mat) - match(selected_name, rownames(mat)) + 1L
  graphics::rect(
    xleft = 0.5,
    ybottom = selected_y - 0.5,
    xright = ncol(mat) + 0.5,
    ytop = selected_y + 0.5,
    border = "black",
    lwd = 2
  )

  cut_positions <- selected$segment_summary$end_index[-nrow(selected$segment_summary)] + 0.5
  if (length(cut_positions) > 0) {
    graphics::abline(v = cut_positions, col = "black", lwd = 1.5, lty = 2)
  }
  summary <- selected$segment_summary
  mids <- (summary$start_index + summary$end_index) / 2
  if (is.null(class_labels)) {
    class_labels <- summary$label
  }
  graphics::text(
    x = mids,
    y = selected_y + 0.42,
    labels = class_labels,
    cex = 0.8,
    font = 2
  )
  invisible(diagnostic)
}

ageTMP_interval_wss_costs <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  sums <- rbind(rep(0, p), apply(x, 2, cumsum))
  sums_sq <- c(0, cumsum(rowSums(x * x)))
  costs <- matrix(Inf, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    j <- i:n
    len <- j - i + 1L
    seg_sums <- sweep(sums[j + 1L, , drop = FALSE], 2, sums[i, ], "-")
    seg_sums_sq <- sums_sq[j + 1L] - sums_sq[i]
    costs[i, j] <- seg_sums_sq - rowSums(seg_sums * seg_sums) / len
  }
  costs
}
