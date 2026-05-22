#' Load manuscript reference-cohort clinical data from STable1
#'
#' @param data_dir Path to the public data directory.
#'
#' @return A data frame containing the `ClinicalTable` sheet from `STable1.xlsx`.
#' @export
ageTMP_load_cdisc_clinical <- function(data_dir = "data") {
  as.data.frame(ageTMP_load_supplement(
    data_dir = data_dir,
    table = "STable1",
    sheet = "ClinicalTable"
  ))
}

#' Load manuscript reference-cohort mutation calls for CPSA modeling
#'
#' The manuscript protein CPSA workflow uses the reference-cohort mutation
#' modeling matrix `cDisc_mutation_10192023.tsv`, which includes both local
#' discovery samples and external reference-labeled samples present in the
#' survival cohort. If that file is unavailable, the function falls back to the
#' `Disc_Mutation` sheet in `STable1.xlsx`, but that reduced supplementary view
#' does not contain the full mutation covariate set used by the manuscript CPSA
#' models.
#'
#' @param data_dir Path to the public data directory.
#' @param genes Character vector of genes to return.
#'
#' @return A numeric sample-by-gene matrix with mutation indicators.
#' @export
ageTMP_load_cdisc_mutation <- function(
  data_dir = "data",
  genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
) {
  aliases <- c(H33A = "H3-3A")
  lookup <- genes
  lookup[lookup %in% names(aliases)] <- aliases[lookup[lookup %in% names(aliases)]]

  mutation_file <- file.path(data_dir, "cDisc_mutation_10192023.tsv")
  if (file.exists(mutation_file)) {
    mutation <- utils::read.delim(mutation_file, sep = "\t", check.names = FALSE)
    gene_col <- "ApprovedGeneSymbol"
  } else {
    mutation <- as.data.frame(ageTMP_load_supplement(
      data_dir = data_dir,
      table = "STable1",
      sheet = "Disc_Mutation"
    ))
    gene_col <- "Gene"
  }

  if (!gene_col %in% names(mutation)) {
    stop("Mutation data are missing required gene column: ", gene_col, call. = FALSE)
  }

  sample_cols <- setdiff(names(mutation), gene_col)
  out <- matrix(NA_real_, nrow = length(sample_cols), ncol = length(genes))
  rownames(out) <- ageTMP_normalize_sample_ids(sample_cols)
  colnames(out) <- paste0(genes, "_mut")

  for (i in seq_along(genes)) {
    row_i <- match(lookup[[i]], mutation[[gene_col]])
    if (is.na(row_i)) {
      next
    }
    calls <- unlist(mutation[row_i, sample_cols], use.names = FALSE)
    numeric_calls <- suppressWarnings(as.numeric(calls))
    if (all(is.na(numeric_calls) & !is.na(calls))) {
      out[, i] <- ifelse(is.na(calls), NA_real_, as.numeric(as.character(calls) != "WT"))
    } else {
      out[, i] <- numeric_calls
    }
  }

  out
}

#' Load manuscript reference-cohort protein data as a gene-by-sample matrix
#'
#' @param data_dir Path to the public data directory.
#' @param collapse Whether to average rows with the same gene symbol.
#'
#' @return A numeric gene-by-sample matrix.
#' @export
ageTMP_load_cdisc_protein_matrix <- function(data_dir = "data", collapse = TRUE) {
  ageTMP_load_feature_matrix(data_dir = data_dir, modality = "protein", collapse = collapse)
}

#' Load a public feature matrix
#'
#' Load one or more public molecular abundance matrices as feature-by-sample
#' numeric matrices. Rows are identified by `ApprovedGeneSymbol` when available,
#' and duplicate feature rows can be averaged to produce one row per molecule.
#'
#' @param data_dir Path to the public data directory.
#' @param modality Molecular modality or modalities accepted by
#'   [ageTMP_load_molecular()]. A single modality returns one matrix; multiple
#'   modalities return a named list of matrices.
#' @param collapse Whether to average rows with the same feature identifier.
#' @param row_id Annotation column to use as the feature identifier.
#'
#' @return A numeric feature-by-sample matrix, or a named list of matrices when
#'   multiple modalities are requested.
#' @export
ageTMP_load_feature_matrix <- function(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho"),
  collapse = TRUE,
  row_id = "ApprovedGeneSymbol"
) {
  allowed <- c("protein", "rna", "glyco", "phospho")
  if (length(modality) > 1) {
    modality <- match.arg(modality, allowed, several.ok = TRUE)
    out <- lapply(modality, function(one_modality) {
      ageTMP_load_feature_matrix(
        data_dir = data_dir,
        modality = one_modality,
        collapse = collapse,
        row_id = row_id
      )
    })
    names(out) <- modality
    return(out)
  }
  modality <- match.arg(modality, allowed)
  raw <- ageTMP_load_molecular(data_dir = data_dir, modality = modality)
  annotation_cols <- switch(
    modality,
    protein = intersect(c("ApprovedGeneSymbol", "Symbol.V5", "OldSymbol", "coding_gene"), names(raw)),
    rna = intersect(c("ApprovedGeneSymbol", "Symbol.V5", "Observation_Rate_HOPE"), names(raw)),
    glyco = seq_len(min(4, ncol(raw))),
    phospho = seq_len(min(9, ncol(raw)))
  )
  if (!row_id %in% names(raw)) {
    row_id <- names(raw)[annotation_cols[1]]
  }
  split <- ageTMP_split_annotation_matrix(raw, annotation_cols = annotation_cols, row_id = row_id)

  if (collapse) {
    ageTMP_collapse_matrix_by_feature(split$matrix, split$annotation[[row_id]])
  } else {
    split$matrix
  }
}

#' Load cDisc RNA matrix
#'
#' @inheritParams ageTMP_load_feature_matrix
#'
#' @return A numeric gene-by-sample RNA matrix.
#' @export
ageTMP_load_cdisc_rna_matrix <- function(data_dir = "data", collapse = TRUE) {
  ageTMP_load_feature_matrix(data_dir = data_dir, modality = "rna", collapse = collapse)
}

#' Load manuscript trajectory QC omit lists
#'
#' The manuscript CPSA workflow uses sex-specific omit lists for features with
#' very flat age-dependent trajectory signal. These omit lists are derived from
#' fitted trajectory matrices by ranking the row-wise standard deviation and
#' using `sd < 0.15` as the default low-dynamic-range threshold.
#'
#' @param rdata_path Path to a `rank_df_rev.RData`-style object.
#' @param type Molecular data type to keep, such as `"protein"`.
#'
#' @return A list with the rank data frame and male/female omit vectors.
#' @export
ageTMP_load_trajectory_qc <- function(rdata_path, type = "protein") {
  if (!file.exists(rdata_path)) {
    stop("Trajectory QC RData file not found: ", rdata_path, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  loaded <- load(rdata_path, envir = env)
  if (!"rank.df" %in% loaded || !"rank.df" %in% ls(env)) {
    stop("`rank.df` was not found in: ", rdata_path, call. = FALSE)
  }

  rank_df <- env$rank.df
  required <- c("gene", "sd", "rank.sd", "sex", "type")
  missing <- setdiff(required, names(rank_df))
  if (length(missing) > 0) {
    stop("Trajectory QC data are missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  rank_df <- rank_df[rank_df$type %in% type, , drop = FALSE]
  male_name <- paste0("m.", type, ".omit")
  female_name <- paste0("f.", type, ".omit")
  list(
    rank_df = rank_df,
    male_omit = if (exists(male_name, envir = env, inherits = FALSE)) get(male_name, envir = env) else character(),
    female_omit = if (exists(female_name, envir = env, inherits = FALSE)) get(female_name, envir = env) else character(),
    source_file = rdata_path,
    type = type
  )
}

#' Load manuscript protein trajectory feature universe
#'
#' @param rdata_path Path to a `protein_rev_tadj62_list.RData`-style object.
#' @param list_index List element containing the manuscript trajectory matrices.
#' @param matrix_name Matrix whose row names define the protein feature universe.
#'
#' @return Character vector of protein features.
#' @export
ageTMP_load_protein_trajectory_features <- function(
  rdata_path,
  list_index = 2,
  matrix_name = "hope.mat0.adj"
) {
  if (!file.exists(rdata_path)) {
    stop("Protein trajectory RData file not found: ", rdata_path, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  load(rdata_path, envir = env)
  if (!exists("m.gene.mat.adj.list", envir = env, inherits = FALSE)) {
    stop("`m.gene.mat.adj.list` was not found in: ", rdata_path, call. = FALSE)
  }
  mat_list <- get("m.gene.mat.adj.list", envir = env)
  if (length(mat_list) < list_index || !matrix_name %in% names(mat_list[[list_index]])) {
    stop("Matrix `", matrix_name, "` was not found at list index ", list_index, ".", call. = FALSE)
  }
  rownames(mat_list[[list_index]][[matrix_name]])
}

#' Prepare reference-cohort clinical covariates for CPSA survival modeling
#'
#' @param clinical Clinical data frame, usually from
#'   [ageTMP_load_reference_clinical()].
#' @param mutation Optional sample-by-gene mutation indicator matrix.
#' @param protein_sample_ids Optional vector of sample IDs available in the protein matrix.
#' @param sex Optional sex filter, `"Male"` or `"Female"`.
#' @param age_classes Age classes to keep.
#' @param max_day Administrative survival truncation in days.
#' @param mutation_na How to handle missing mutation covariates after matching
#'   mutation calls to clinical samples. Use `"zero"` to treat missing calls as
#'   wild type/absent, or `"preserve"` to keep missing values as `NA`. The
#'   manuscript protein CPSA workflow used sex-specific behavior: missing
#'   mutation covariates were zero-filled in the male run and preserved in the
#'   female run.
#' @param require_complete_mutation For public reference clinical tables that
#'   contain `has.complete.mut`, restrict to samples with complete mutation
#'   covariates. This matches the manuscript reference-cohort CPSA model frame.
#'
#' @return A data frame with survival outcome and model covariates.
#' @export
ageTMP_prepare_cdisc_cpsa_survival <- function(
  clinical,
  mutation = NULL,
  protein_sample_ids = NULL,
  sex = NULL,
  age_classes = c("PED", "ADO", "YA", "ADULT"),
  max_day = 2000,
  mutation_na = c("zero", "preserve"),
  require_complete_mutation = TRUE
) {
  mutation_na <- match.arg(mutation_na)
  required <- c(
    "id", "cDisc_os", "cDisc_os_status", "cDisc_pfs", "cDisc_pfs_status",
    "cDisc_clinical_status_at_collection_event", "cDisc_treat_status",
    "cDisc_First_Diagnosis", "cDisc_age", "cDisc_Gender",
    "cDisc_WHO_Grade", "cDisc_tumor_loc"
  )
  missing <- setdiff(required, names(clinical))
  if (length(missing) > 0) {
    stop("Clinical data are missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  id <- ageTMP_normalize_sample_ids(clinical$id)
  os <- suppressWarnings(as.numeric(clinical$cDisc_os))
  pfs <- suppressWarnings(as.numeric(clinical$cDisc_pfs))
  os_status <- suppressWarnings(as.numeric(clinical$cDisc_os_status))
  pfs_status <- suppressWarnings(as.numeric(clinical$cDisc_pfs_status))

  os_status <- ifelse(os > max_day, 0, os_status)
  pfs_status <- ifelse(pfs > max_day, 0, pfs_status)
  os_surv <- ifelse(os > max_day, max_day, os)
  pfs_surv <- ifelse(pfs > max_day, max_day, pfs)

  out <- data.frame(
    id = id,
    days = os_surv,
    pfs.days = pfs_surv,
    os.status = os_status,
    pfs.status = pfs_status,
    dasc = clinical$cDisc_clinical_status_at_collection_event,
    treat_status = clinical$cDisc_treat_status,
    First.Diagnosis = clinical$cDisc_First_Diagnosis,
    age = suppressWarnings(as.numeric(clinical$cDisc_age)),
    Gender = clinical$cDisc_Gender,
    Grade = as.numeric(suppressWarnings(as.numeric(clinical$cDisc_WHO_Grade)) %in% 3),
    Midline = clinical$cDisc_tumor_loc %in% "Midline",
    Cortical = clinical$cDisc_tumor_loc %in% "Cortical",
    stringsAsFactors = FALSE
  )

  out$age_class_new <- factor(
    cut(out$age, breaks = c(0, 15, 26, 40, 62), include.lowest = TRUE,
        labels = c("PED", "ADO", "YA", "ADULT")),
    levels = c("PED", "ADO", "YA", "ADULT")
  )

  if (!is.null(mutation)) {
    mutation <- as.matrix(mutation)
    mutation <- mutation[match(out$id, rownames(mutation)), , drop = FALSE]
    if (identical(mutation_na, "zero")) {
      mutation[is.na(mutation)] <- 0
    }
    out <- cbind(out, as.data.frame(mutation, check.names = FALSE))
  }

  if (!is.null(protein_sample_ids)) {
    out <- out[out$id %in% ageTMP_normalize_sample_ids(protein_sample_ids), , drop = FALSE]
  }
  if (isTRUE(require_complete_mutation) && "has.complete.mut" %in% names(clinical)) {
    out <- out[clinical$has.complete.mut[match(out$id, ageTMP_normalize_sample_ids(clinical$id))] %in% TRUE, , drop = FALSE]
  }
  if (!is.null(sex)) {
    out <- out[out$Gender %in% sex, , drop = FALSE]
  }

  out <- out[out$age_class_new %in% age_classes, , drop = FALSE]
  out <- out[!out$dasc %in% "Deceased-due to disease", , drop = FALSE]
  out <- out[out$treat_status %in% "Treatment naive", , drop = FALSE]
  out <- out[out$First.Diagnosis %in% "Yes", , drop = FALSE]
  out <- out[!is.na(out$days) & !is.na(out$os.status), , drop = FALSE]
  rownames(out) <- out$id
  out
}

#' Create a CPSA model specification
#'
#' CPSA is feature-matrix agnostic: each row of `feature_matrix` is a molecular
#' feature or score, and each column is a sample. This specification object makes
#' the survival outcome, age-stratum interaction, and clinical/molecular
#' adjustment model explicit so users can adapt CPSA to different applications.
#'
#' @param time_col Column in `survival_data` containing survival time.
#' @param event_col Column in `survival_data` containing event status.
#' @param feature_col Temporary column name used for the current feature.
#' @param age_class_col Column containing age-class labels. Set to `NULL` with
#'   `age_reference = NULL` and `strata = NULL` to fit a standard feature Cox
#'   model without feature-by-age-class interactions.
#' @param age_reference Reference level for `age_class_col`.
#' @param strata Age strata to test by dropping feature terms.
#' @param include_all_combined_test Whether to add an omnibus `ALL.comb.*`
#'   likelihood-ratio test that drops the main feature and all feature-by-strata
#'   interaction terms. Legacy two-level glyco ADO scripts did not create this
#'   column.
#' @param base_covariates Clinical covariates included in every model.
#' @param adjustment_covariates Additional covariates, such as mutation calls.
#' @param fdr_method Multiple-testing adjustment method passed to
#'   [stats::p.adjust()].
#' @param remove_sparse_covariates Whether to drop binary covariates with very
#'   few positive values in the current model frame.
#' @param sparse_covariates Optional character vector limiting which covariates
#'   are eligible for sparse binary checks. If `NULL`, only covariates with no
#'   more than two observed values are checked.
#' @param sparse_positive_level Positive value used for sparse binary checks.
#' @param sparse_min_count Drop sparse covariates with counts less than or equal
#'   to this value.
#' @param sparse_count_method How to count sparse positive values. `"numeric"`
#'   coerces binary/logical covariates to numeric before counting. `"factor"`
#'   preserves the legacy manuscript behavior of counting
#'   `table(factor(x, levels = c(0, 1)))[2]`, which drops logical covariates
#'   such as `Cortical` and `Midline` in the glyco ADO scripts.
#' @param scale_features Whether to row-center and row-scale features inside
#'   the model engine. Set to `FALSE` when the supplied matrix has already been
#'   scaled in a manuscript-specific way.
#' @param engine Cox model fitting engine. `"coxph"` uses the formula-based
#'   [survival::coxph()] path and is the default for manuscript replication.
#'   `"coxph.fit"` uses the same model matrix with [survival::coxph.fit()] to
#'   reduce formula overhead in larger feature scans.
#'
#' @return A list describing the CPSA Cox model.
#' @export
ageTMP_cpsa_spec <- function(
  time_col = "days",
  event_col = "os.status",
  feature_col = "feature_score",
  age_class_col = "age_class_new",
  age_reference = "ADULT",
  strata = c("PED", "ADO", "YA"),
  include_all_combined_test = TRUE,
  base_covariates = c("age", "Grade", "Cortical", "Midline"),
  adjustment_covariates = c("IDH1_mut", "TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut"),
  fdr_method = "BY",
  remove_sparse_covariates = TRUE,
  sparse_covariates = NULL,
  sparse_positive_level = 1,
  sparse_min_count = 3,
  sparse_count_method = c("numeric", "factor"),
  scale_features = TRUE,
  engine = c("coxph", "coxph.fit")
) {
  sparse_count_method <- match.arg(sparse_count_method)
  engine <- match.arg(engine)
  list(
    time_col = time_col,
    event_col = event_col,
    feature_col = feature_col,
    age_class_col = age_class_col,
    age_reference = age_reference,
    strata = strata,
    include_all_combined_test = include_all_combined_test,
    base_covariates = base_covariates,
    adjustment_covariates = adjustment_covariates,
    fdr_method = fdr_method,
    remove_sparse_covariates = remove_sparse_covariates,
    sparse_covariates = sparse_covariates,
    sparse_positive_level = sparse_positive_level,
    sparse_min_count = sparse_min_count,
    sparse_count_method = sparse_count_method,
    scale_features = scale_features,
    engine = engine
  )
}

#' Fit reference-cohort CPSA Cox models
#'
#' This is the modality-neutral CPSA model engine. It accepts any numeric
#' feature-by-sample matrix, including protein abundance, RNA expression,
#' glycopeptide abundance, phosphosite abundance, pathway scores, or other
#' molecular feature scores.
#'
#' @param feature_matrix Numeric feature-by-sample matrix.
#' @param survival_data Prepared survival covariates.
#' @param features Optional subset of features to fit.
#' @param spec CPSA model specification from [ageTMP_cpsa_spec()].
#' @param trajectory_sd Optional trajectory dynamic-range object used to screen
#'   features before CPSA fitting. Accepted inputs are a named numeric vector of
#'   per-feature trajectory standard deviations, a feature-by-age/sample matrix
#'   whose row-wise standard deviation should be used, or a data frame such as
#'   the output of [ageTMP_rank_trajectory_sd()].
#' @param trajectory_sd_min Minimum trajectory standard deviation required to
#'   keep a feature. Use `NULL` to disable trajectory-SD screening.
#' @param trajectory_sd_feature_col Feature column in a `trajectory_sd` data
#'   frame.
#' @param trajectory_sd_value_col Standard-deviation column in a `trajectory_sd`
#'   data frame.
#' @param trajectory_sd_group_cols Optional grouping columns in a `trajectory_sd`
#'   data frame. When supplied, `trajectory_sd_keep` controls whether a feature
#'   must pass in all groups or any group.
#' @param trajectory_sd_keep Whether grouped trajectory-SD screening keeps
#'   features passing in `"any_group"` or `"all_groups"`.
#' @param n_cores Number of parallel worker processes on non-Windows systems.
#' @param progress Whether to print fitting progress.
#' @param progress_every Number of features per progress chunk.
#'
#' @return A data frame with coefficients, p-values, adjusted FDR values, and
#'   signed age-stratum statistics.
#' @export
ageTMP_fit_reference_cpsa <- function(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  spec = ageTMP_cpsa_spec(),
  trajectory_sd = NULL,
  trajectory_sd_min = NULL,
  trajectory_sd_feature_col = "feature",
  trajectory_sd_value_col = "sd",
  trajectory_sd_group_cols = NULL,
  trajectory_sd_keep = c("any_group", "all_groups"),
  n_cores = 1,
  progress = FALSE,
  progress_every = 500
) {
  trajectory_sd_keep <- match.arg(trajectory_sd_keep)
  feature_matrix <- as.matrix(feature_matrix)
  colnames(feature_matrix) <- ageTMP_normalize_sample_ids(colnames(feature_matrix))
  features <- intersect(features, rownames(feature_matrix))
  features <- ageTMP_screen_cpsa_features_by_trajectory_sd(
    features = features,
    trajectory_sd = trajectory_sd,
    sd_min = trajectory_sd_min,
    feature_col = trajectory_sd_feature_col,
    value_col = trajectory_sd_value_col,
    group_cols = trajectory_sd_group_cols,
    keep = trajectory_sd_keep
  )
  if (length(features) == 0) {
    stop("No features remain after applying the CPSA feature filters.", call. = FALSE)
  }
  sample_ids <- intersect(colnames(feature_matrix), rownames(survival_data))
  if (length(sample_ids) < 5) {
    stop("Fewer than five shared feature/survival samples are available.", call. = FALSE)
  }

  scaled <- feature_matrix[features, sample_ids, drop = FALSE]
  if (isTRUE(spec$scale_features)) {
    scaled <- t(scale(t(scaled)))
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1) {
    n_cores <- 1
  }
  progress_every <- as.integer(progress_every)
  if (is.na(progress_every) || progress_every < 1) {
    progress_every <- 500
  }

  meta <- survival_data[sample_ids, , drop = FALSE]
  fit_one <- function(i) {
    ageTMP_fit_one_cpsa_feature(
      score = as.numeric(scaled[i, ]),
      meta = meta,
      spec = spec
    )
  }

  indices <- seq_len(nrow(scaled))
  chunks <- split(indices, ceiling(seq_along(indices) / progress_every))
  fit_list <- vector("list", length(indices))
  if (isTRUE(progress)) {
    message("ageTMP CPSA fit 0/", length(indices))
  }
  for (chunk in chunks) {
    chunk_fits <- if (n_cores > 1 && .Platform$OS.type != "windows") {
      parallel::mclapply(chunk, fit_one, mc.cores = min(n_cores, length(chunk)))
    } else {
      lapply(chunk, fit_one)
    }
    fit_list[chunk] <- chunk_fits
    if (isTRUE(progress)) {
      message("ageTMP CPSA fit ", max(chunk), "/", length(indices))
    }
  }

  ageTMP_collect_cpsa_fits(rownames(scaled), fit_list, spec = spec)
}

ageTMP_screen_cpsa_features_by_trajectory_sd <- function(
  features,
  trajectory_sd = NULL,
  sd_min = NULL,
  feature_col = "feature",
  value_col = "sd",
  group_cols = NULL,
  keep = c("any_group", "all_groups")
) {
  keep <- match.arg(keep)
  if (is.null(sd_min)) {
    return(features)
  }
  sd_min <- suppressWarnings(as.numeric(sd_min))
  if (length(sd_min) != 1 || is.na(sd_min)) {
    stop("`trajectory_sd_min` must be a single numeric value or `NULL`.", call. = FALSE)
  }
  if (is.null(trajectory_sd)) {
    stop("`trajectory_sd` must be supplied when `trajectory_sd_min` is used.", call. = FALSE)
  }

  if (is.matrix(trajectory_sd) || is.data.frame(trajectory_sd) && all(vapply(trajectory_sd, is.numeric, logical(1)))) {
    mat <- as.matrix(trajectory_sd)
    if (is.null(rownames(mat))) {
      stop("Matrix-like `trajectory_sd` input must have feature row names.", call. = FALSE)
    }
    sd_values <- stats::setNames(apply(mat, 1, stats::sd, na.rm = TRUE), rownames(mat))
    return(intersect(features, names(sd_values)[!is.na(sd_values) & sd_values >= sd_min]))
  }

  if (is.atomic(trajectory_sd) && is.numeric(trajectory_sd)) {
    if (is.null(names(trajectory_sd))) {
      stop("Numeric `trajectory_sd` input must be named by feature.", call. = FALSE)
    }
    return(intersect(features, names(trajectory_sd)[!is.na(trajectory_sd) & trajectory_sd >= sd_min]))
  }

  trajectory_sd <- data.frame(trajectory_sd, stringsAsFactors = FALSE)
  required <- c(feature_col, value_col, group_cols)
  missing <- setdiff(required, names(trajectory_sd))
  if (length(missing) > 0) {
    stop("`trajectory_sd` is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  sd_rank <- trajectory_sd
  names(sd_rank)[names(sd_rank) == feature_col] <- "feature"
  names(sd_rank)[names(sd_rank) == value_col] <- "sd"
  if (length(group_cols) > 0) {
    kept <- ageTMP_filter_trajectory_sd(
      sd_rank = sd_rank,
      sd_min = sd_min,
      group_cols = group_cols,
      keep = keep
    )
  } else {
    sd_rank$passes_sd <- !is.na(sd_rank$sd) & sd_rank$sd >= sd_min
    pass_by_feature <- tapply(sd_rank$passes_sd, sd_rank$feature, any)
    kept <- names(pass_by_feature)[pass_by_feature]
  }
  intersect(features, kept)
}

#' Fit reference-cohort protein CPSA Cox models
#'
#' This fits the reference-cohort CPSA model used for manuscript protein
#' features:
#' a Cox model with scaled protein abundance, protein-by-age-class interactions,
#' age, tumor covariates, and mutation covariates. Combined age-stratum p-values
#' are obtained by dropping the stratum-specific protein term from the full model.
#'
#' @param protein_matrix Numeric gene-by-sample protein matrix.
#' @param survival_data Prepared survival covariates from
#'   [ageTMP_prepare_reference_cpsa_survival()].
#' @param genes Optional subset of genes to fit.
#' @param age_reference Reference age class for the interaction model.
#' @param mutation_genes Mutation gene names without `_mut`.
#'
#' @return A data frame with reference-cohort CPSA statistics.
#' @export
ageTMP_fit_cdisc_protein_cpsa <- function(
  protein_matrix,
  survival_data,
  genes = rownames(protein_matrix),
  age_reference = "ADULT",
  mutation_genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
) {
  ageTMP_fit_reference_cpsa(
    feature_matrix = protein_matrix,
    survival_data = survival_data,
    features = genes,
    spec = ageTMP_cpsa_spec(
      feature_col = "score2",
      age_reference = age_reference,
      adjustment_covariates = paste0(mutation_genes, "_mut")
    )
  )
}

ageTMP_fit_one_cpsa_feature <- function(
  score,
  meta,
  spec = ageTMP_cpsa_spec()
) {
  feature_col <- spec$feature_col
  age_class_mode <- !is.null(spec$age_class_col) &&
    !is.null(spec$age_reference) &&
    length(spec$strata) > 0
  data <- data.frame(meta, check.names = FALSE)
  data[[feature_col]] <- score
  data <- data[!is.na(data[[spec$time_col]]) & !is.na(data[[spec$event_col]]), , drop = FALSE]
  if ("First.Diagnosis" %in% names(data)) {
    data <- data[data$First.Diagnosis %in% "Yes", , drop = FALSE]
  }

  adjustment_covariates <- intersect(spec$adjustment_covariates, names(data))
  base_covariates <- intersect(spec$base_covariates, names(data))

  keep <- unique(c(
    feature_col, if (age_class_mode) spec$age_class_col, spec$time_col, spec$event_col,
    base_covariates, adjustment_covariates
  ))
  data <- data[, keep, drop = FALSE]
  if (age_class_mode) {
    data[[spec$age_class_col]] <- stats::relevel(
      factor(data[[spec$age_class_col]], levels = unique(c(spec$strata, spec$age_reference))),
      ref = spec$age_reference
    )
  }
  data$SurvObj <- survival::Surv(data[[spec$time_col]], data[[spec$event_col]])

  if (isTRUE(spec$remove_sparse_covariates)) {
    sparse_check <- function(x) {
      if (identical(spec$sparse_count_method, "factor")) {
        counts <- table(factor(x, levels = c(0, 1)))
        return(as.numeric(counts[["1"]]) <= spec$sparse_min_count)
      }
      numeric_x <- suppressWarnings(as.numeric(x))
      sum(numeric_x %in% spec$sparse_positive_level, na.rm = TRUE) <= spec$sparse_min_count
    }
    sparse_covariates <- intersect(c(base_covariates, adjustment_covariates), names(data))
    if (!is.null(spec$sparse_covariates)) {
      sparse_covariates <- intersect(sparse_covariates, spec$sparse_covariates)
    } else {
      sparse_covariates <- sparse_covariates[
        vapply(data[, sparse_covariates, drop = FALSE], function(x) {
          length(unique(stats::na.omit(x))) <= 2
        }, logical(1))
      ]
    }
    remove <- sparse_covariates[vapply(data[, sparse_covariates, drop = FALSE], sparse_check, logical(1))]
    if (length(remove) > 0) {
      data <- data[, !names(data) %in% remove, drop = FALSE]
    }
  }

  model_terms <- c(
    feature_col,
    if (age_class_mode) paste0(spec$age_class_col, ":", feature_col),
    intersect(base_covariates, names(data)),
    intersect(adjustment_covariates, names(data))
  )
  formula <- stats::reformulate(model_terms, response = "SurvObj")

  old_na <- getOption("na.action")
  on.exit(options(na.action = old_na), add = TRUE)
  options(na.action = "na.pass")
  design <- data.frame(stats::model.matrix(formula, data, na.action = "na.pass")[, -1, drop = FALSE])
  options(na.action = "na.omit")
  design$SurvObj <- data$SurvObj

  full_terms <- setdiff(names(design), "SurvObj")
  fit_full <- ageTMP_safe_cox_model(full_terms, design, full = TRUE, engine = spec$engine)
  if (is.null(fit_full)) {
    return(list(anova = NULL, coef = NULL))
  }

  pvals <- numeric(0)
  if (age_class_mode) {
    age_terms <- c(
      stats::setNames(
        paste0(feature_col, ".", spec$age_class_col, spec$strata),
        spec$strata
      )
    )
    reduced_terms <- stats::setNames(
      lapply(spec$strata, function(stratum) c(feature_col, age_terms[[stratum]])),
      spec$strata
    )
    if (isTRUE(spec$include_all_combined_test)) {
      reduced_terms$ALL <- c(feature_col, unname(age_terms))
    }

    pvals <- rep(NA_real_, length(reduced_terms))
    names(pvals) <- paste0(names(reduced_terms), ".comb.p")
    for (nm in names(reduced_terms)) {
      keep_terms <- setdiff(setdiff(names(design), "SurvObj"), reduced_terms[[nm]])
      reduced <- design[, c(keep_terms, "SurvObj"), drop = FALSE]
      fit_reduced <- ageTMP_safe_cox_model(keep_terms, reduced, full = FALSE, engine = spec$engine)
      if (!is.null(fit_reduced)) {
        pvals[[paste0(nm, ".comb.p")]] <- ageTMP_cox_lrt_p(fit_full, fit_reduced)
      }
    }
  }

  list(anova = pvals, coef = ageTMP_cox_coef_table(fit_full))
}

ageTMP_safe_cox_model <- function(terms, data, full = FALSE, engine = c("coxph", "coxph.fit")) {
  engine <- match.arg(engine)
  if (identical(engine, "coxph")) {
    return(ageTMP_safe_coxph(stats::reformulate(terms, response = "SurvObj"), data, full = full))
  }
  ageTMP_safe_coxph_fit(terms, data, full = full)
}

ageTMP_safe_coxph <- function(formula, data, full = FALSE) {
  tryCatch(
    survival::coxph(
      formula,
      data = data,
      control = ageTMP_coxph_control(full = full)
    ),
    warning = function(w) NULL,
    error = function(e) NULL
  )
}

ageTMP_coxph_control <- function(full = FALSE) {
  survival::coxph.control(
    eps = 1e-09,
    toler.chol = .Machine$double.eps ^ 0.75,
    iter.max = if (full) 1000000 else 1000,
    toler.inf = sqrt(1e-09),
    outer.max = if (full) 100 else 10,
    timefix = TRUE
  )
}

ageTMP_safe_coxph_fit <- function(terms, data, full = FALSE) {
  if (length(terms) == 0) {
    return(NULL)
  }
  x <- as.matrix(data[, terms, drop = FALSE])
  keep <- stats::complete.cases(x) & stats::complete.cases(data$SurvObj)
  x <- x[keep, , drop = FALSE]
  y <- data$SurvObj[keep]
  if (nrow(x) < 2) {
    return(NULL)
  }
  storage.mode(x) <- "double"
  fit <- tryCatch(
    survival::coxph.fit(
      x = x,
      y = y,
      strata = rep(0L, nrow(x)),
      offset = rep(0, nrow(x)),
      init = NULL,
      control = ageTMP_coxph_control(full = full),
      weights = rep(1, nrow(x)),
      method = "efron",
      rownames = rownames(data)[keep],
      resid = FALSE,
      nocenter = NULL
    ),
    warning = function(w) NULL,
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$coefficients) || is.null(fit$var)) {
    return(NULL)
  }
  names(fit$coefficients) <- colnames(x)
  dimnames(fit$var) <- list(colnames(x), colnames(x))
  fit$terms <- terms
  class(fit) <- c("ageTMP_coxph_fit", class(fit))
  fit
}

ageTMP_cox_coef_table <- function(fit) {
  if (inherits(fit, "ageTMP_coxph_fit")) {
    coef <- fit$coefficients
    se <- sqrt(diag(fit$var))
    z <- coef / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    out <- cbind(
      coef = coef,
      `exp(coef)` = exp(coef),
      `se(coef)` = se,
      z = z,
      `Pr(>|z|)` = p
    )
    rownames(out) <- names(coef)
    return(out)
  }
  summary(fit)$coef
}

ageTMP_cox_lrt_p <- function(fit_full, fit_reduced) {
  if (inherits(fit_full, "ageTMP_coxph_fit") && inherits(fit_reduced, "ageTMP_coxph_fit")) {
    df <- sum(!is.na(fit_full$coefficients)) - sum(!is.na(fit_reduced$coefficients))
    if (df <= 0) {
      return(NA_real_)
    }
    chisq <- 2 * (fit_full$loglik[2] - fit_reduced$loglik[2])
    return(stats::pchisq(chisq, df = df, lower.tail = FALSE))
  }
  tryCatch(
    as.matrix(stats::anova(fit_full, fit_reduced))[2, 4],
    error = function(e) NA_real_
  )
}

ageTMP_collect_cpsa_fits <- function(genes, fit_list, spec = ageTMP_cpsa_spec()) {
  set_matrix_colnames <- function(x, value, label) {
    if (ncol(x) != length(value)) {
      stop(
        label, " has ", ncol(x), " column(s) but ", length(value), " column name(s).",
        call. = FALSE
      )
    }
    colnames(x) <- value
    x
  }
  paste0_or_empty <- function(x, suffix) {
    if (length(x) == 0) character() else paste0(x, suffix)
  }

  ok <- vapply(fit_list, function(x) !is.null(x$coef) && !is.null(x$anova), logical(1))
  genes <- genes[ok]
  fit_list <- fit_list[ok]
  if (length(fit_list) == 0) {
    return(data.frame(pathway = character(), check.names = FALSE))
  }

  coef_names <- unique(unlist(lapply(fit_list, function(x) rownames(x$coef))))
  coef_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(coef_names))
  p_mat <- se_mat <- z_mat <- coef_mat
  rownames(coef_mat) <- rownames(p_mat) <- rownames(se_mat) <- rownames(z_mat) <- genes
  coef_mat <- set_matrix_colnames(coef_mat, paste0(coef_names, ".coef"), "coef_mat")
  se_mat <- set_matrix_colnames(se_mat, paste0(coef_names, ".se"), "se_mat")
  z_mat <- set_matrix_colnames(z_mat, paste0(coef_names, ".zscore"), "z_mat")
  p_mat <- set_matrix_colnames(p_mat, paste0(coef_names, ".p"), "p_mat")

  age_class_mode <- !is.null(spec$age_class_col) &&
    !is.null(spec$age_reference) &&
    length(spec$strata) > 0
  comb_names <- character()
  if (age_class_mode) {
    comb_names <- spec$strata
    if (isTRUE(spec$include_all_combined_test)) {
      comb_names <- c(comb_names, "ALL")
    }
  }
  comb_p <- matrix(NA_real_, nrow = length(genes), ncol = length(comb_names))
  comb_p <- set_matrix_colnames(comb_p, paste0_or_empty(comb_names, ".comb.p"), "comb_p")
  rownames(comb_p) <- genes

  for (i in seq_along(fit_list)) {
    co <- fit_list[[i]]$coef
    rn <- rownames(co)
    coef_mat[i, paste0(rn, ".coef")] <- co[, 1]
    se_mat[i, paste0(rn, ".se")] <- co[, 3]
    z_mat[i, paste0(rn, ".zscore")] <- co[, 4]
    p_mat[i, paste0(rn, ".p")] <- co[, 5]
    comb_p[i, names(fit_list[[i]]$anova)] <- fit_list[[i]]$anova
  }

  fdr_mat <- p_mat
  for (j in seq_len(ncol(p_mat))) {
    fdr_mat[, j] <- stats::p.adjust(p_mat[, j], method = spec$fdr_method)
  }
  rownames(fdr_mat) <- rownames(p_mat)
  fdr_mat <- set_matrix_colnames(fdr_mat, sub("\\.p$", ".fdr", colnames(p_mat)), "fdr_mat")
  comb_fdr <- comb_p
  for (j in seq_len(ncol(comb_p))) {
    comb_fdr[, j] <- stats::p.adjust(comb_p[, j], method = spec$fdr_method)
  }
  rownames(comb_fdr) <- rownames(comb_p)
  comb_fdr <- set_matrix_colnames(comb_fdr, sub("\\.p$", ".fdr", colnames(comb_p)), "comb_fdr")

  score_coef <- coef_mat[, paste0(spec$feature_col, ".coef")]
  if (age_class_mode) {
    int_coef <- function(term) {
      col <- paste0(spec$feature_col, ".", spec$age_class_col, term, ".coef")
      if (col %in% colnames(coef_mat)) coef_mat[, col] else rep(0, nrow(coef_mat))
    }
    comb_coef <- do.call(
      cbind,
      stats::setNames(
        lapply(spec$strata, function(stratum) score_coef + int_coef(stratum)),
        paste0(spec$strata, ".comb.coef")
      )
    )
    if (isTRUE(spec$include_all_combined_test)) {
      comb_coef <- cbind(
        comb_coef,
        ALL.comb.coef = score_coef + rowSums(do.call(cbind, lapply(spec$strata, int_coef)))
      )
    }
    comb_coef <- matrix(
      as.numeric(comb_coef),
      nrow = nrow(comb_p),
      ncol = ncol(comb_p),
      dimnames = list(rownames(comb_p), sub("\\.p$", ".coef", colnames(comb_p)))
    )
  } else {
    comb_coef <- matrix(NA_real_, nrow = length(genes), ncol = 0)
    rownames(comb_coef) <- genes
  }

  signed_p <- p_mat
  signed_p[,] <- log10(p_mat) * -1 * sign(coef_mat)
  signed_p <- set_matrix_colnames(signed_p, sub("\\.p$", ".signed.pval", colnames(p_mat)), "signed_p")
  signed_fdr <- fdr_mat
  signed_fdr[,] <- log10(fdr_mat) * -1 * sign(coef_mat)
  signed_fdr <- set_matrix_colnames(signed_fdr, sub("\\.fdr$", ".signed.fdr", colnames(fdr_mat)), "signed_fdr")
  comb_signed_p <- comb_p
  if (ncol(comb_p) > 0) {
    comb_signed_p[,] <- log10(comb_p) * -1 * sign(comb_coef)
  }
  comb_signed_p <- set_matrix_colnames(comb_signed_p, sub("\\.p$", ".signed.p", colnames(comb_p)), "comb_signed_p")
  comb_signed_fdr <- comb_fdr
  if (ncol(comb_fdr) > 0) {
    comb_signed_fdr[,] <- log10(comb_fdr) * -1 * sign(comb_coef)
  }
  comb_signed_fdr <- set_matrix_colnames(comb_signed_fdr, sub("\\.fdr$", ".signed.fdr", colnames(comb_fdr)), "comb_signed_fdr")

  data.frame(
    pathway = genes,
    signed_fdr,
    signed_p,
    coef_mat,
    comb_coef,
    p_mat,
    fdr_mat,
    z_mat,
    se_mat,
    comb_p,
    comb_fdr,
    comb_signed_p,
    comb_signed_fdr,
    check.names = FALSE
  )
}

#' Build manuscript reference-cohort protein CPSA columns for STable4 comparison
#'
#' @param data_dir Path to the public data directory.
#' @param genes Optional gene subset for testing.
#' @param mode Reproducibility mode. `"manuscript"` preserves the original
#'   sex-specific mutation missingness behavior used by the manuscript protein
#'   CPSA scripts: male mutation NAs are set to zero and female mutation NAs are
#'   preserved. `"standardized"` zero-fills missing mutation calls for both
#'   sexes.
#'
#' @return A data frame with manuscript reference-cohort signed log10p columns
#'   for male and female. Published `STable4` column names retain the `cdisc`
#'   label because that is the manuscript source-table convention.
#' @export
ageTMP_build_sa_protein_cdisc <- function(
  data_dir = "data",
  genes = NULL,
  mode = c("manuscript", "standardized")
) {
  mode <- match.arg(mode)
  clinical <- ageTMP_load_cdisc_clinical(data_dir)
  mutation <- ageTMP_load_cdisc_mutation(data_dir)
  protein <- ageTMP_load_cdisc_protein_matrix(data_dir)
  if (!is.null(genes)) {
    protein <- protein[intersect(genes, rownames(protein)), , drop = FALSE]
  }

  fit_one_sex <- function(sex) {
    mutation_na <- if (identical(mode, "manuscript") && identical(sex, "Female")) {
      "preserve"
    } else {
      "zero"
    }
    surv <- ageTMP_prepare_cdisc_cpsa_survival(
      clinical = clinical,
      mutation = mutation,
      protein_sample_ids = colnames(protein),
      sex = sex,
      mutation_na = mutation_na
    )
    fit <- ageTMP_fit_cdisc_protein_cpsa(protein, surv)
    keep <- c("pathway", "PED.comb.signed.p", "ADO.comb.signed.p", "YA.comb.signed.p")
    fit[, keep, drop = FALSE]
  }

  male <- fit_one_sex("Male")
  female <- fit_one_sex("Female")
  out <- merge(male, female, by = "pathway", all = TRUE, suffixes = c(".male", ".female"))
  names(out)[names(out) == "pathway"] <- "gene"
  names(out) <- sub("comb\\.signed\\.p", "comb.signed.log10p.cdisc", names(out))
  out[order(out$gene), , drop = FALSE]
}

#' Build manuscript STable4 discovery-cohort CPSA columns
#'
#' This modality-neutral helper fits the cDisc/discovery CPSA model to a public
#' feature matrix and returns the signed log10 p-value columns used in STable4.
#'
#' @param data_dir Path to the public data directory.
#' @param modality Molecular modality accepted by [ageTMP_load_feature_matrix()].
#' @param features Optional feature subset.
#' @param mode Reproducibility mode. `"manuscript"` preserves the sex-specific
#'   mutation missingness behavior used in the manuscript CPSA scripts.
#'
#' @return A data frame with signed log10 p-value columns for male and female.
#' @export
ageTMP_build_sa_discovery_cpsa <- function(
  data_dir = "data",
  modality = c("protein", "rna", "glyco", "phospho"),
  features = NULL,
  mode = c("manuscript", "standardized")
) {
  modality <- match.arg(modality)
  mode <- match.arg(mode)
  clinical <- ageTMP_load_cdisc_clinical(data_dir)
  mutation <- ageTMP_load_cdisc_mutation(data_dir)
  feature_matrix <- ageTMP_load_feature_matrix(data_dir, modality = modality)
  if (!is.null(features)) {
    feature_matrix <- feature_matrix[intersect(features, rownames(feature_matrix)), , drop = FALSE]
  }

  fit_one_sex <- function(sex) {
    mutation_na <- if (identical(mode, "manuscript") && identical(sex, "Female")) {
      "preserve"
    } else {
      "zero"
    }
    surv <- ageTMP_prepare_cdisc_cpsa_survival(
      clinical = clinical,
      mutation = mutation,
      protein_sample_ids = colnames(feature_matrix),
      sex = sex,
      mutation_na = mutation_na
    )
    fit <- ageTMP_fit_reference_cpsa(
      feature_matrix = feature_matrix,
      survival_data = surv,
      spec = ageTMP_cpsa_spec(feature_col = "score2")
    )
    keep <- c("pathway", "PED.comb.signed.p", "ADO.comb.signed.p", "YA.comb.signed.p")
    fit[, keep, drop = FALSE]
  }

  male <- fit_one_sex("Male")
  female <- fit_one_sex("Female")
  out <- merge(male, female, by = "pathway", all = TRUE, suffixes = c(".male", ".female"))
  names(out)[names(out) == "pathway"] <- "gene"
  names(out) <- sub("comb\\.signed\\.p", "comb.signed.log10p.cdisc", names(out))
  out[order(out$gene), , drop = FALSE]
}

#' Load reference-cohort clinical data
#'
#' Public-facing alias for [ageTMP_load_cdisc_clinical()]. In the manuscript
#' source files, the external/reference survival cohort is labeled `cDisc`; in
#' the package API, "reference cohort" is the preferred general term.
#'
#' @inheritParams ageTMP_load_cdisc_clinical
#'
#' @return A clinical data frame.
#' @export
ageTMP_load_reference_clinical <- function(data_dir = "data") {
  as.data.frame(ageTMP_load_supplement(
    data_dir = data_dir,
    table = "STable1",
    sheet = "Ref_ClinicalTable"
  ))
}

#' Load discovery-cohort clinical data
#'
#' Public-facing alias for [ageTMP_load_cdisc_clinical()].
#'
#' @inheritParams ageTMP_load_cdisc_clinical
#'
#' @return A clinical data frame.
#' @export
ageTMP_load_discovery_clinical <- function(data_dir = "data") {
  ageTMP_load_cdisc_clinical(data_dir = data_dir)
}

#' Load reference-cohort mutation covariates
#'
#' Public-facing alias for [ageTMP_load_cdisc_mutation()]. The returned matrix is
#' suitable for reference-cohort CPSA adjustment covariates.
#'
#' @inheritParams ageTMP_load_cdisc_mutation
#'
#' @return A numeric sample-by-gene mutation indicator matrix.
#' @export
ageTMP_load_reference_mutation <- function(
  data_dir = "data",
  genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
) {
  ageTMP_load_cdisc_mutation(data_dir = data_dir, genes = genes)
}

#' Load reference-cohort protein matrix
#'
#' Public-facing alias for [ageTMP_load_cdisc_protein_matrix()]. The manuscript
#' implementation currently supports protein CPSA first; the naming leaves room
#' for RNA, glyco, phosphosite, and other feature matrices to use the same
#' reference-cohort vocabulary.
#'
#' @inheritParams ageTMP_load_cdisc_protein_matrix
#'
#' @return A numeric gene-by-sample protein matrix.
#' @export
ageTMP_load_reference_protein_matrix <- function(data_dir = "data", collapse = TRUE) {
  ageTMP_load_cdisc_protein_matrix(data_dir = data_dir, collapse = collapse)
}

#' Prepare a reference-cohort CPSA survival model frame
#'
#' Public-facing alias for [ageTMP_prepare_cdisc_cpsa_survival()]. It constructs
#' the survival outcome and adjustment covariates used by the protein CPSA
#' reference-cohort model.
#'
#' @inheritParams ageTMP_prepare_cdisc_cpsa_survival
#'
#' @return A data frame with survival outcome and model covariates.
#' @export
ageTMP_prepare_reference_cpsa_survival <- function(
  clinical,
  mutation = NULL,
  protein_sample_ids = NULL,
  sex = NULL,
  age_classes = c("PED", "ADO", "YA", "ADULT"),
  max_day = 2000,
  mutation_na = c("zero", "preserve"),
  require_complete_mutation = TRUE
) {
  mutation_na <- match.arg(mutation_na)
  clinical <- data.frame(clinical, stringsAsFactors = FALSE)

  if (all(c("cDisc_os", "cDisc_os_status") %in% names(clinical))) {
    return(ageTMP_prepare_cdisc_cpsa_survival(
      clinical = clinical,
      mutation = mutation,
      protein_sample_ids = protein_sample_ids,
      sex = sex,
      age_classes = age_classes,
      max_day = max_day,
      mutation_na = mutation_na,
      require_complete_mutation = FALSE
    ))
  }

  required <- c("id", "days", "os.status", "age", "age_class_name", "Gender", "Diag", "First.Diagnosis")
  missing <- setdiff(required, names(clinical))
  if (length(missing) > 0) {
    stop("Reference clinical data are missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  days <- suppressWarnings(as.numeric(clinical$days))
  status <- suppressWarnings(as.numeric(clinical$os.status))
  status <- ifelse(days > max_day, 0, status)
  days <- ifelse(days > max_day, max_day, days)

  out <- data.frame(
    id = ageTMP_normalize_sample_ids(clinical$id),
    days = days,
    os.status = status,
    age = suppressWarnings(as.numeric(clinical$age)),
    age_class_new = clinical$age_class_name,
    Gender = clinical$Gender,
    First.Diagnosis = clinical$First.Diagnosis,
    gbm = as.numeric(clinical$Diag %in% "Glioblastoma"),
    mdg = as.numeric(clinical$Diag %in% "Midline"),
    stringsAsFactors = FALSE
  )

  mutation_cols <- grep("_mut$", names(clinical), value = TRUE)
  if (length(mutation_cols) > 0) {
    mutation_df <- clinical[, mutation_cols, drop = FALSE]
    mutation_df <- as.data.frame(lapply(mutation_df, function(x) suppressWarnings(as.numeric(x))))
    if (identical(mutation_na, "zero")) {
      mutation_df[is.na(mutation_df)] <- 0
    }
    out <- cbind(out, mutation_df)
  }

  if (!is.null(mutation)) {
    mutation <- as.matrix(mutation)
    mutation <- mutation[match(out$id, rownames(mutation)), , drop = FALSE]
    if (identical(mutation_na, "zero")) {
      mutation[is.na(mutation)] <- 0
    }
    out <- cbind(out, as.data.frame(mutation, check.names = FALSE))
  }

  if (!is.null(protein_sample_ids)) {
    out <- out[out$id %in% ageTMP_normalize_sample_ids(protein_sample_ids), , drop = FALSE]
  }
  if (isTRUE(require_complete_mutation) && "has.complete.mut" %in% names(clinical)) {
    complete_mut <- clinical$has.complete.mut[match(out$id, ageTMP_normalize_sample_ids(clinical$id))]
    out <- out[complete_mut %in% TRUE, , drop = FALSE]
  }
  if (!is.null(sex)) {
    out <- out[out$Gender %in% sex, , drop = FALSE]
  }

  out <- out[out$age_class_new %in% age_classes, , drop = FALSE]
  out <- out[out$First.Diagnosis %in% "Yes", , drop = FALSE]
  out <- out[!is.na(out$days) & !is.na(out$os.status), , drop = FALSE]
  rownames(out) <- out$id
  out
}

#' Fit reference-cohort protein CPSA models
#'
#' Public-facing alias for [ageTMP_fit_cdisc_protein_cpsa()]. This is the model
#' engine for reference-cohort protein CPSA.
#'
#' @inheritParams ageTMP_fit_cdisc_protein_cpsa
#'
#' @return A data frame with coefficients, p-values, BY-adjusted FDR values, and
#'   signed combined age-stratum statistics.
#' @export
ageTMP_fit_reference_protein_cpsa <- function(
  protein_matrix,
  survival_data,
  genes = rownames(protein_matrix),
  age_reference = "ADULT",
  mutation_genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
) {
  ageTMP_fit_cdisc_protein_cpsa(
    protein_matrix = protein_matrix,
    survival_data = survival_data,
    genes = genes,
    age_reference = age_reference,
    mutation_genes = mutation_genes
  )
}

#' Build manuscript STable4 protein reference-cohort CPSA columns
#'
#' Public-facing alias for [ageTMP_build_sa_protein_cdisc()]. It returns the
#' manuscript `SA-Protein-cDisc-Ref` reference-cohort signed log10 p-value
#' columns, preserving the published column naming convention.
#'
#' @inheritParams ageTMP_build_sa_protein_cdisc
#'
#' @return A data frame with manuscript reference-cohort signed log10 p-value
#'   columns for male and female.
#' @export
ageTMP_build_sa_protein_reference <- function(
  data_dir = "data",
  genes = NULL,
  mode = c("manuscript", "standardized")
) {
  ageTMP_build_sa_protein_cdisc(data_dir = data_dir, genes = genes, mode = mode)
}
