#' Rank trajectory features by standard deviation
#'
#' Compute the within-group standard deviation of fitted trajectory values.
#' This makes explicit the manuscript QC idea behind `rank_df_rev.RData`: some
#' downstream analyses exclude features with very flat age-dependent tumor
#' trajectories, commonly using `sd < 0.15` as the omit threshold.
#'
#' @param trajectory Long-format trajectory data frame.
#' @param feature_col Column containing feature names.
#' @param value_col Column containing trajectory values, usually `fit`.
#' @param group_cols Columns defining independent rankings. For manuscript
#'   protein DTT-style filtering this is typically sex, and optionally tissue
#'   when normal and tumor trajectories are present in the same table.
#' @param tissue_col Optional tissue column used with `tissue`.
#' @param tissue Optional tissue value to retain before ranking, such as
#'   `"Tumor"`. Use `NULL` to rank all rows supplied.
#'
#' @return A data frame with grouping columns, `feature`, `sd`, and `rank_sd`.
#' @export
ageTMP_rank_trajectory_sd <- function(
  trajectory,
  feature_col = "feature",
  value_col = "fit",
  group_cols = c("sex"),
  tissue_col = "tissue",
  tissue = NULL
) {
  trajectory <- data.frame(trajectory, stringsAsFactors = FALSE)
  required <- unique(c(feature_col, value_col, group_cols))
  if (!is.null(tissue)) {
    required <- unique(c(required, tissue_col))
  }
  missing <- setdiff(required, names(trajectory))
  if (length(missing) > 0) {
    stop(
      "`trajectory` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(tissue)) {
    trajectory <- trajectory[trajectory[[tissue_col]] %in% tissue, , drop = FALSE]
  }
  if (nrow(trajectory) == 0) {
    stop("No trajectory rows remain after filtering.", call. = FALSE)
  }

  split_key <- if (length(group_cols) == 0) {
    factor(rep("all", nrow(trajectory)))
  } else {
    interaction(trajectory[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  }

  out <- lapply(split(trajectory, split_key), function(df) {
    values <- split(df[[value_col]], df[[feature_col]])
    sd_values <- vapply(values, stats::sd, numeric(1), na.rm = TRUE)
    rank_values <- rank(sd_values, ties.method = "average", na.last = "keep")
    group_data <- if (length(group_cols) == 0) {
      data.frame(stringsAsFactors = FALSE)
    } else {
      unique(df[, group_cols, drop = FALSE])[1, , drop = FALSE]
    }
    data.frame(
      group_data[rep(1, length(sd_values)), , drop = FALSE],
      feature = names(sd_values),
      sd = as.numeric(sd_values),
      rank_sd = as.numeric(rank_values),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Filter trajectory features by standard deviation
#'
#' Apply an explicit dynamic-range filter to the output of
#' [ageTMP_rank_trajectory_sd()]. The default keeps features with
#' `sd >= 0.15`, matching the omit-list threshold used in the manuscript
#' `rank_df_rev.RData` workflow.
#'
#' @param sd_rank Data frame returned by [ageTMP_rank_trajectory_sd()].
#' @param sd_min Minimum standard deviation required to keep a feature.
#' @param group_cols Columns defining groups that must pass the threshold.
#' @param keep One of `"all_groups"`, `"any_group"`, or `"per_group"`.
#'   `"all_groups"` keeps features passing in every requested group;
#'   `"any_group"` keeps features passing in at least one group; `"per_group"`
#'   returns the group-specific rows that pass.
#'
#' @return A character vector of kept features for `"all_groups"` or
#'   `"any_group"`, and a filtered data frame for `"per_group"`.
#' @export
ageTMP_filter_trajectory_sd <- function(
  sd_rank,
  sd_min = 0.15,
  group_cols = c("sex"),
  keep = c("all_groups", "any_group", "per_group")
) {
  keep <- match.arg(keep)
  sd_rank <- data.frame(sd_rank, stringsAsFactors = FALSE)
  required <- unique(c("feature", "sd", group_cols))
  missing <- setdiff(required, names(sd_rank))
  if (length(missing) > 0) {
    stop(
      "`sd_rank` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  sd_rank$passes_sd <- !is.na(sd_rank$sd) & sd_rank$sd >= sd_min
  if (keep == "per_group") {
    return(sd_rank[sd_rank$passes_sd, , drop = FALSE])
  }

  pass_by_feature <- tapply(sd_rank$passes_sd, sd_rank$feature, if (keep == "all_groups") all else any)
  names(pass_by_feature)[pass_by_feature]
}
