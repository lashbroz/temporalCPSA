#' Select dynamic TMP features
#'
#' Select rows from a temporal molecular profile matrix by row-wise standard
#' deviation, then optionally row-center and row-scale the selected matrix.
#'
#' @param tmp_matrix Numeric feature-by-sample TMP matrix.
#' @param proportion Proportion of rows to keep after ranking by standard
#'   deviation.
#' @param max_features Maximum number of rows to keep.
#' @param scale_rows Whether to row-center and row-scale the selected matrix.
#'
#' @return A numeric matrix containing selected dynamic TMP rows.
#' @export
ageTMP_select_dynamic_tmp_features <- function(
  tmp_matrix,
  proportion = 0.5,
  max_features = 5000,
  scale_rows = TRUE
) {
  tmp_matrix <- as.matrix(tmp_matrix)
  if (is.null(rownames(tmp_matrix)) || is.null(colnames(tmp_matrix))) {
    stop("`tmp_matrix` must have row and column names.", call. = FALSE)
  }
  keep_row <- rowSums(!is.na(tmp_matrix)) >= 2
  tmp_matrix <- tmp_matrix[keep_row, , drop = FALSE]
  if (nrow(tmp_matrix) == 0) {
    stop("No TMP rows have at least two non-missing values.", call. = FALSE)
  }

  row_sd <- apply(tmp_matrix, 1, stats::sd, na.rm = TRUE)
  row_sd[is.na(row_sd)] <- -Inf
  proportion <- suppressWarnings(as.numeric(proportion))
  max_features <- suppressWarnings(as.integer(max_features))
  if (length(proportion) != 1 || is.na(proportion) || proportion <= 0 || proportion > 1) {
    stop("`proportion` must be a single number in (0, 1].", call. = FALSE)
  }
  if (length(max_features) != 1 || is.na(max_features) || max_features < 1) {
    stop("`max_features` must be a positive integer.", call. = FALSE)
  }

  n_keep <- min(max_features, max(1L, round(nrow(tmp_matrix) * proportion)))
  selected <- tmp_matrix[order(row_sd, decreasing = TRUE)[seq_len(n_keep)], , drop = FALSE]
  if (isTRUE(scale_rows)) {
    selected <- t(scale(t(selected)))
    selected <- selected[rowSums(is.na(selected)) == 0, , drop = FALSE]
  }
  selected
}

#' Combine TMP matrices for multi-omic sample clustering
#'
#' Select dynamic rows within each molecular data type, restrict to common
#' samples, and stack the matrices for sample-level clustering. Row names are
#' prefixed by modality so that repeated feature identifiers remain unique.
#'
#' @param tmp_matrices Named list of numeric feature-by-sample TMP matrices.
#' @inheritParams ageTMP_select_dynamic_tmp_features
#' @param common_samples Optional sample IDs to require. If `NULL`, the
#'   intersection of matrix column names is used.
#'
#' @return A numeric stacked TMP matrix.
#' @export
ageTMP_combine_tmp_matrices <- function(
  tmp_matrices,
  proportion = 0.5,
  max_features = 5000,
  common_samples = NULL,
  scale_rows = TRUE
) {
  if (!is.list(tmp_matrices) || length(tmp_matrices) == 0) {
    stop("`tmp_matrices` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(tmp_matrices)) || any(!nzchar(names(tmp_matrices)))) {
    names(tmp_matrices) <- paste0("modality", seq_along(tmp_matrices))
  }
  tmp_matrices <- lapply(tmp_matrices, as.matrix)
  sample_sets <- lapply(tmp_matrices, colnames)
  if (is.null(common_samples)) {
    common_samples <- Reduce(intersect, sample_sets)
  } else {
    common_samples <- Reduce(intersect, c(list(ageTMP_normalize_sample_ids(common_samples)), sample_sets))
  }
  if (length(common_samples) < 2) {
    stop("Fewer than two common samples are available across TMP matrices.", call. = FALSE)
  }

  selected <- lapply(names(tmp_matrices), function(modality) {
    mat <- tmp_matrices[[modality]][, common_samples, drop = FALSE]
    mat <- ageTMP_select_dynamic_tmp_features(
      mat,
      proportion = proportion,
      max_features = max_features,
      scale_rows = scale_rows
    )
    mat <- as.matrix(mat)
    if (nrow(mat) == 0) {
      stop("No dynamic TMP rows remain for modality `", modality, "`.", call. = FALSE)
    }
    rn <- rownames(mat)
    if (is.null(rn) || length(rn) != nrow(mat)) {
      rn <- paste0("feature", seq_len(nrow(mat)))
    }
    dimnames(mat) <- list(paste(modality, make.unique(as.character(rn)), sep = "::"), colnames(mat))
    mat
  })
  names(selected) <- names(tmp_matrices)
  do.call(rbind, selected)
}

#' Cluster samples from temporal molecular profiles
#'
#' Cluster samples using the stacked TMP matrix produced by
#' [ageTMP_combine_tmp_matrices()]. When `ConsensusClusterPlus` is available,
#' this uses the consensus k-means workflow used by the manuscript clustering
#' scripts. Otherwise, it falls back to deterministic k-means for lightweight
#' development and testing.
#'
#' @param tmp_matrix Numeric feature-by-sample TMP matrix.
#' @param k_max Maximum number of clusters to evaluate.
#' @param reps Number of consensus resampling repetitions.
#' @param p_item Sample resampling proportion for consensus clustering.
#' @param p_feature Feature resampling proportion for consensus clustering.
#' @param algorithm Clustering algorithm passed to `ConsensusClusterPlus`.
#' @param distance Distance metric passed to `ConsensusClusterPlus`.
#' @param seed Random seed.
#' @param title Output title/directory prefix used by `ConsensusClusterPlus`.
#' @param plot Plot option passed to `ConsensusClusterPlus`.
#' @param use_consensus Whether to use `ConsensusClusterPlus` when available.
#'
#' @return A list with the input `clustme`, cluster assignments by k, and the
#'   raw consensus result when applicable.
#' @export
ageTMP_cluster_tmp_samples <- function(
  tmp_matrix,
  k_max = 10,
  reps = 300,
  p_item = 1,
  p_feature = 1,
  algorithm = "km",
  distance = "euclidean",
  seed = 1234,
  title = "ageTMP_consensus_tmp",
  plot = NULL,
  use_consensus = TRUE
) {
  tmp_matrix <- as.matrix(tmp_matrix)
  if (is.null(rownames(tmp_matrix)) || is.null(colnames(tmp_matrix))) {
    stop("`tmp_matrix` must have row and column names.", call. = FALSE)
  }
  tmp_matrix <- tmp_matrix[rowSums(is.na(tmp_matrix)) == 0, colSums(is.na(tmp_matrix)) == 0, drop = FALSE]
  if (nrow(tmp_matrix) < 2 || ncol(tmp_matrix) < 2) {
    stop("`tmp_matrix` must contain at least two complete rows and columns.", call. = FALSE)
  }
  k_max <- as.integer(k_max)
  if (is.na(k_max) || k_max < 2) {
    stop("`k_max` must be at least 2.", call. = FALSE)
  }

  if (isTRUE(use_consensus) && requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    raw <- ConsensusClusterPlus::ConsensusClusterPlus(
      tmp_matrix,
      maxK = k_max,
      reps = reps,
      pItem = p_item,
      pFeature = p_feature,
      title = title,
      clusterAlg = algorithm,
      distance = distance,
      seed = seed,
      plot = plot
    )
    clusters <- lapply(raw[2:length(raw)], function(x) x$consensusClass)
    names(clusters) <- paste0("k", seq_along(clusters) + 1L)
    return(list(clustme = tmp_matrix, clusters = clusters, raw = raw))
  }

  set.seed(seed)
  clusters <- lapply(seq.int(2, k_max), function(k) {
    stats::kmeans(t(tmp_matrix), centers = k, nstart = 50)$cluster
  })
  names(clusters) <- paste0("k", seq.int(2, k_max))
  list(clustme = tmp_matrix, clusters = clusters, raw = NULL)
}

#' Summarize age cutpoints from TMP sample clusters
#'
#' Order clustered samples by age, relabel clusters by first appearance along
#' the age axis, and summarize each cluster's age range. This mirrors the final
#' interpretation step in the manuscript age-class clustering workflow.
#'
#' @param clusters Named cluster-assignment vector.
#' @param clinical Clinical metadata.
#' @param sample_col Clinical sample ID column.
#' @param age_col Clinical age column.
#'
#' @return A list containing relabeled clusters and an age-range summary.
#' @export
ageTMP_summarize_tmp_age_clusters <- function(
  clusters,
  clinical,
  sample_col = "id",
  age_col = "age"
) {
  clinical <- data.frame(clinical, stringsAsFactors = FALSE)
  required <- c(sample_col, age_col)
  missing <- setdiff(required, names(clinical))
  if (length(missing) > 0) {
    stop("`clinical` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ids <- ageTMP_normalize_sample_ids(clinical[[sample_col]])
  ages <- suppressWarnings(as.numeric(clinical[[age_col]]))
  names(ages) <- ids
  clusters <- clusters[names(clusters) %in% names(ages)]
  clusters <- clusters[!is.na(ages[names(clusters)])]
  clusters <- clusters[order(ages[names(clusters)])]
  if (length(clusters) == 0) {
    stop("No clustered samples match clinical ages.", call. = FALSE)
  }
  relabeled <- factor(clusters, levels = unique(clusters), labels = seq_along(unique(clusters)))
  names(relabeled) <- names(clusters)
  age_range <- stats::aggregate(
    age ~ cluster,
    data = data.frame(cluster = relabeled, age = ages[names(relabeled)]),
    FUN = function(x) paste(range(x, na.rm = TRUE), collapse = "-")
  )
  names(age_range)[names(age_range) == "age"] <- "age_range"
  list(clusters = relabeled, age_range = age_range)
}
