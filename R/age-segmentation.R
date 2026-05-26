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
#' @return A list with `matrix`, `ordered_samples`, per-`k` `fits`, and a
#'   `summary` table containing the segmentation score and mean silhouette width
#'   for each candidate `k`.
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
  tmp_for_silhouette <- as.matrix(tmp_matrix)
  if (isTRUE(scale_rows)) {
    tmp_for_silhouette <- t(scale(t(tmp_for_silhouette)))
  }
  tmp_for_silhouette <- tmp_for_silhouette[rowSums(is.na(tmp_for_silhouette)) == 0, , drop = FALSE]
  colnames(tmp_for_silhouette) <- ageTMP_normalize_sample_ids(colnames(tmp_for_silhouette))

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
  tmp_for_silhouette <- tmp_for_silhouette[, ordered_ids, drop = FALSE]
  sample_dist <- stats::dist(t(tmp_for_silhouette))
  mean_silhouette <- vapply(fits, function(fit) {
    if (!requireNamespace("cluster", quietly = TRUE) || fit$k < 2) {
      return(NA_real_)
    }
    classes <- as.integer(fit$classes[ordered_ids])
    if (length(unique(classes)) < 2) {
      return(NA_real_)
    }
    mean(cluster::silhouette(classes, sample_dist)[, "sil_width"], na.rm = TRUE)
  }, numeric(1))
  summary <- data.frame(
    k = k_values,
    score = vapply(fits, `[[`, numeric(1), "score"),
    mean_silhouette = mean_silhouette,
    stringsAsFactors = FALSE
  )
  list(matrix = mat, ordered_samples = fits[[1]]$ordered_samples, fits = fits, summary = summary)
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
#' @param selected_k Which fitted `k` to overlay. If `NULL`, defaults to the
#'   fitted `k` with the largest mean silhouette when available, otherwise the
#'   largest fitted `k`.
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
#' @param show_split_hierarchy Whether to annotate the dominant contiguous
#'   split and next supported splits from the fitted segmentation summary.
#' @param split_label_k Optional integer k values used for hierarchy labels.
#'   By default, the selected k is labeled as dominant and the next best
#'   non-selected k by mean silhouette is labeled as the next supported split.
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
  show_split_hierarchy = TRUE,
  split_label_k = NULL,
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
    if ("summary" %in% names(diagnostic) && "mean_silhouette" %in% names(diagnostic$summary) &&
        any(!is.na(diagnostic$summary$mean_silhouette))) {
      selected_name <- paste0("k", diagnostic$summary$k[which.max(diagnostic$summary$mean_silhouette)])
    } else {
      selected_name <- utils::tail(names(diagnostic$fits), 1)
    }
  } else {
    selected_name <- paste0("k", as.integer(selected_k))
  }
  if (!selected_name %in% names(diagnostic$fits)) {
    stop("Selected `k` was not fit in `diagnostic`.", call. = FALSE)
  }
  selected <- diagnostic$fits[[selected_name]]
  age <- diagnostic$ordered_samples$age

  old_par <- graphics::par(c("mar", "xpd"))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(if (isTRUE(show_age_axis)) 5 else 2, 5, 4, 6), xpd = NA)
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
    ticks <- seq(floor(min(age, na.rm = TRUE)), ceiling(max(age, na.rm = TRUE)), by = 1)
    tick_pos <- stats::approx(age, seq_along(age), xout = ticks, rule = 2, ties = "ordered")$y
    keep <- ticks >= min(age, na.rm = TRUE) & ticks <= max(age, na.rm = TRUE)
    graphics::axis(side = 1, at = tick_pos[keep], labels = ticks[keep], las = 2, cex.axis = 0.6)
  }

  if ("summary" %in% names(diagnostic) && "mean_silhouette" %in% names(diagnostic$summary)) {
    diagnostic_summary <- diagnostic$summary
    label_by_k <- stats::setNames(
      ifelse(
        is.na(diagnostic_summary$mean_silhouette),
        "NA",
        sprintf("%.2f", diagnostic_summary$mean_silhouette)
      ),
      paste0("k", diagnostic_summary$k)
    )
    row_labels <- label_by_k[rownames(mat)]
    graphics::mtext("Mean\nsil.", side = 4, line = 2.7, at = nrow(mat) + 0.4, cex = 0.75, font = 2)
    graphics::text(
      x = ncol(mat) + 3.5,
      y = seq_len(nrow(mat)),
      labels = rev(row_labels),
      cex = 0.75,
      adj = 0
    )
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
      graphics::segments(
        x0 = suggested_pos,
        y0 = 0.5,
        x1 = suggested_pos,
        y1 = nrow(mat) + 0.5,
        col = "#111111",
        lwd = 2,
        lty = 1
      )
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
    graphics::segments(
      x0 = cut_positions,
      y0 = 0.5,
      x1 = cut_positions,
      y1 = nrow(mat) + 0.5,
      col = "black",
      lwd = 2.4,
      lty = 2
    )
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

  if (isTRUE(show_split_hierarchy) && "summary" %in% names(diagnostic)) {
    summary_df <- diagnostic$summary
    if (is.null(split_label_k)) {
      split_label_k <- as.integer(sub("^k", "", selected_name))
      if ("mean_silhouette" %in% names(summary_df)) {
        remaining <- summary_df[summary_df$k != split_label_k & !is.na(summary_df$mean_silhouette), , drop = FALSE]
        if (nrow(remaining) > 0) {
          split_label_k <- c(split_label_k, remaining$k[which.max(remaining$mean_silhouette)])
        }
      }
    }
    split_label_k <- unique(split_label_k[!is.na(split_label_k)])
    split_labels <- c("Dominant split", "Next supported split", "Additional supported split")
    label_rows <- seq_len(min(length(split_label_k), length(split_labels)))
    for (i in label_rows) {
      fit_name <- paste0("k", split_label_k[[i]])
      if (!fit_name %in% names(diagnostic$fits)) {
        next
      }
      fit <- diagnostic$fits[[fit_name]]
      if (nrow(fit$segment_summary) < 2) {
        next
      }
      pos <- fit$segment_summary$end_index[-nrow(fit$segment_summary)] + 0.5
      y_pos <- nrow(mat) - match(fit_name, rownames(mat)) + 1L
      graphics::segments(
        x0 = pos,
        y0 = y_pos - 0.42,
        x1 = pos,
        y1 = y_pos + 0.42,
        col = "#111111",
        lwd = if (i == 1L) 3.2 else 2.4
      )
      label <- split_labels[[i]]
      if (i == 1L && length(pos) == 1L) {
        label <- paste0(label, "\n~", round(age[min(length(age), max(1, floor(pos)))], 0), " yr")
      }
      graphics::text(
        x = min(max(pos), ncol(mat) + 1.2),
        y = y_pos + 0.58,
        labels = label,
        cex = if (i == 1L) 0.78 else 0.68,
        font = 2,
        adj = c(0, 0),
        col = "#111111"
      )
    }
  }
  invisible(diagnostic)
}

#' Plot age-segmentation depth from contiguous AD-TMP structure
#'
#' Draw a simplified depth plot for contiguous age-segmentation diagnostics.
#' Each row shows the age cutpoints selected for a candidate `k`, with darker
#' bars and thicker cutpoint markers emphasizing deeper supported splits. This
#' plot is intended for interpretation and communication of the segmentation
#' hierarchy, while [ageTMP_plot_segment_diagnostic()] provides the fuller
#' ordered diagnostic heatmap.
#'
#' @param diagnostic Output of [ageTMP_segment_diagnostic_matrix()].
#' @param xlim Numeric length-two age range for the display.
#' @param main Plot title.
#' @param dominant_label Label for the selected strongest split.
#' @param support_labels Optional named character vector mapping k values to
#'   support-level labels.
#' @param show_sample_rug Whether to show sample ages as a rug track.
#' @param show_legend Whether to show a compact explanatory legend.
#'
#' @return Invisibly returns `diagnostic`.
#' @export
ageTMP_plot_segment_depth <- function(
  diagnostic,
  xlim = c(0, 80),
  main = "Data-Driven Age-Segmentation Depth from AD-TMP Structure",
  dominant_label = "max silhouette: strongest split",
  support_labels = NULL,
  show_sample_rug = TRUE,
  show_legend = TRUE
) {
  if (!is.list(diagnostic) || !all(c("fits", "summary", "ordered_samples") %in% names(diagnostic))) {
    stop("`diagnostic` must be the output of `ageTMP_segment_diagnostic_matrix()`.", call. = FALSE)
  }
  summary_df <- diagnostic$summary
  if (!all(c("k", "mean_silhouette") %in% names(summary_df))) {
    stop("`diagnostic$summary` must contain `k` and `mean_silhouette`.", call. = FALSE)
  }
  cutpoints_by_k <- lapply(summary_df$k, function(k) {
    fit <- diagnostic$fits[[paste0("k", k)]]
    if (is.null(fit) || nrow(fit$segment_summary) < 2) {
      return(numeric())
    }
    fit$cutpoints
  })
  if (is.null(support_labels)) {
    support_labels <- c(
      `2` = "deepest",
      `3` = "next broad",
      `4` = "finer",
      `5` = "intermediate",
      `6` = "intermediate",
      `7` = "finer",
      `8` = "finer"
    )
  }
  k_vals <- summary_df$k
  y <- rev(seq_along(k_vals)) * 1.18
  max_y <- max(y)
  selected_i <- if ("selected_by_max_mean_silhouette" %in% names(summary_df) &&
      any(isTRUE(summary_df$selected_by_max_mean_silhouette))) {
    which(summary_df$selected_by_max_mean_silhouette)[[1]]
  } else {
    which.max(summary_df$mean_silhouette)
  }
  dominant_values <- cutpoints_by_k[[selected_i]]
  dominant <- if (summary_df$k[[selected_i]] == 2L && length(dominant_values) == 1L) {
    dominant_values[[1]]
  } else {
    NA_real_
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(5.8, 4.8, 3.5, 4.8), xpd = FALSE)
  graphics::plot(
    NA,
    xlim = c(xlim[[1]] - 18, xlim[[2]] + 16),
    ylim = c(-0.75, max_y + 1.0),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    cex.main = 1.15,
    font.main = 2
  )
  graphics::axis(1, at = seq(xlim[[1]], xlim[[2]], by = 5), cex.axis = 0.85)
  graphics::mtext("Age in years", side = 1, line = 3.4, cex = 0.9)
  graphics::text(xlim[[1]] - 16.0, max_y + 0.62, labels = "support level", adj = 0, cex = 0.78, font = 2)

  green_dark <- "#006d2c"
  green_light <- "#d9f0d3"
  label_col <- c(deepest = "black", `next broad` = "#238b45", finer = "#555555", intermediate = "#555555")
  for (i in seq_along(k_vals)) {
    k <- k_vals[[i]]
    yi <- y[[i]]
    cuts <- cutpoints_by_k[[i]]
    cuts <- cuts[!is.na(cuts)]
    boundaries <- c(xlim[[1]], cuts, xlim[[2]])
    n_seg <- length(boundaries) - 1L
    seg_cols <- grDevices::colorRampPalette(c(green_light, green_dark))(max(2L, n_seg))
    for (j in seq_len(n_seg)) {
      graphics::segments(
        boundaries[j],
        yi,
        boundaries[j + 1L],
        yi,
        col = seg_cols[j],
        lwd = 10,
        lend = "butt"
      )
    }
    if (length(cuts) > 0) {
      is_dominant_cut <- !is.na(dominant) & abs(cuts - dominant) < 0.1
      lwd <- ifelse(is_dominant_cut, 5, 3)
      graphics::segments(cuts, yi - 0.16, cuts, yi + 0.16, col = "#222222", lwd = lwd, lend = "round")
      cut_labels <- round(cuts, 1)
      cut_labels[i == selected_i & is_dominant_cut] <- ""
      graphics::text(cuts, yi + 0.58, labels = cut_labels, srt = 55, cex = 0.58, col = "#444444")
    }
    level <- unname(support_labels[as.character(k)])
    if (is.na(level)) level <- "candidate"
    level_color <- label_col[[level]]
    if (is.null(level_color)) level_color <- "#555555"
    graphics::text(xlim[[1]] - 16.0, yi, labels = level, adj = 0, cex = 0.72, col = level_color)
    graphics::text(xlim[[1]] - 0.2, yi, labels = paste0("K=", k), adj = 1, cex = 0.78)
    graphics::text(
      xlim[[2]] + 4,
      yi,
      labels = sprintf("mean sil. %.2f", summary_df$mean_silhouette[[i]]),
      adj = 0,
      cex = 0.78,
      col = if (i == selected_i) "black" else if (k == 3) "#238b45" else "#333333"
    )
  }

  if (!is.na(dominant)) {
    dom_y <- y[[selected_i]]
    graphics::segments(dominant, dom_y - 0.22, dominant, dom_y + 0.24, col = "black", lwd = 7, lend = "round")
    dom_age <- round(dominant)
    graphics::text(dominant + 0.15, dom_y + 0.66, labels = dom_age, srt = 55, font = 2, cex = 0.75)
    graphics::arrows(dominant + 5.3, dom_y + 0.62, dominant + 0.25, dom_y + 0.25, length = 0.08, lwd = 1.5)
    graphics::text(
      dominant + 5.6,
      dom_y + 0.76,
      labels = dominant_label,
      adj = 0,
      cex = 0.82,
      font = 2
    )
  }

  if (isTRUE(show_sample_rug) && "age" %in% names(diagnostic$ordered_samples)) {
    ages <- sort(unique(round(as.numeric(diagnostic$ordered_samples$age), 3)))
    rug_y <- 0.62
    graphics::segments(ages, rug_y - 0.05, ages, rug_y + 0.05, col = "#c9c9c9", lwd = 0.8)
    graphics::text(xlim[[1]], rug_y - 0.26, labels = "sample ages", adj = 0, cex = 0.65, col = "#555555")
  }

  if (isTRUE(show_legend)) {
    legend_y <- -0.32
    legend_x <- c(xlim[[1]] - 6, xlim[[1]] + 27, xlim[[1]] + 58)
    graphics::segments(legend_x[[1]], legend_y, legend_x[[1]] + 3, legend_y, col = green_dark, lwd = 5, lend = "round")
    graphics::text(legend_x[[1]] + 4.4, legend_y, labels = "darker = deeper supported split", adj = 0, cex = 0.62)
    if (!is.na(dominant)) {
      graphics::segments(legend_x[[2]], legend_y, legend_x[[2]] + 3, legend_y, col = "black", lwd = 7, lend = "round")
      graphics::text(legend_x[[2]] + 4.4, legend_y, labels = paste0("dominant ~", round(dominant), " cutpoint"), adj = 0, cex = 0.62)
      graphics::segments(legend_x[[3]], legend_y, legend_x[[3]] + 3, legend_y, col = "#222222", lwd = 3, lend = "round")
      graphics::text(legend_x[[3]] + 4.4, legend_y, labels = "additional cutpoints", adj = 0, cex = 0.62)
    } else {
      graphics::segments(legend_x[[2]], legend_y, legend_x[[2]] + 3, legend_y, col = "#222222", lwd = 3, lend = "round")
      graphics::text(legend_x[[2]] + 4.4, legend_y, labels = "candidate cutpoints", adj = 0, cex = 0.62)
    }
  }
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
