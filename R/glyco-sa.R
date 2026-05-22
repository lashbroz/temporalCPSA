#' Load public cDisc glycopeptide data
#'
#' Unlike gene-level protein/RNA loaders, glycopeptide survival analyses operate
#' at the glycopeptide row level. This loader preserves `Gene.Sequence` row
#' identifiers and returns the full public glycopeptide annotation.
#'
#' @param data_dir Path to the public data directory.
#'
#' @return A list with `matrix` and `annotation` elements.
#' @export
ageTMP_load_cdisc_glyco_matrix <- function(data_dir = "data") {
  raw <- ageTMP_load_molecular(data_dir = data_dir, modality = "glyco")
  annotation_cols <- intersect(
    c(
      "Modification", "ApprovedGeneSymbol", "Gene", "Sequence",
      "Gene.Sequence", "GlycanType", "gene.type"
    ),
    names(raw)
  )
  required <- c("Modification", "Gene", "Gene.Sequence")
  missing <- setdiff(required, annotation_cols)
  if (length(missing) > 0) {
    stop("Glyco source is missing annotation column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  split <- ageTMP_split_annotation_matrix(
    raw,
    annotation_cols = annotation_cols,
    row_id = "Modification"
  )

  rownames(split$matrix) <- split$annotation$Gene.Sequence
  rownames(split$annotation) <- split$annotation$Gene.Sequence
  split
}

#' Build source-derived protein-adjusted glycopeptide matrix
#'
#' This reproduces the manuscript adjusted glycopeptide matrix by residualizing
#' each glycopeptide abundance vector against the matched gene-level protein
#' abundance vector over shared cDisc samples.
#'
#' @param data_dir Path to the public data directory.
#' @param glyco Optional output from [ageTMP_load_cdisc_glyco_matrix()].
#' @param protein_matrix Optional gene-by-sample protein matrix.
#'
#' @return A list with `matrix` and `annotation` elements.
#' @export
ageTMP_build_cdisc_glyco_adjusted_matrix <- function(
  data_dir = "data",
  glyco = NULL,
  protein_matrix = NULL
) {
  if (is.null(glyco)) {
    glyco <- ageTMP_load_cdisc_glyco_matrix(data_dir = data_dir)
  }
  if (is.null(protein_matrix)) {
    protein_matrix <- ageTMP_load_cdisc_protein_matrix(data_dir = data_dir, collapse = TRUE)
  }

  glyco_matrix <- as.matrix(glyco$matrix)
  glyco_annotation <- data.frame(glyco$annotation, stringsAsFactors = FALSE)
  protein_matrix <- as.matrix(protein_matrix)
  colnames(glyco_matrix) <- ageTMP_normalize_sample_ids(colnames(glyco_matrix))
  colnames(protein_matrix) <- ageTMP_normalize_sample_ids(colnames(protein_matrix))

  sample_ids <- intersect(colnames(glyco_matrix), colnames(protein_matrix))
  if (length(sample_ids) < 3) {
    stop("Fewer than three shared glyco/protein samples are available.", call. = FALSE)
  }

  protein_rows <- match(glyco_annotation$Gene, rownames(protein_matrix))
  adjusted <- vector("list", nrow(glyco_matrix))
  names(adjusted) <- rownames(glyco_matrix)

  for (i in seq_len(nrow(glyco_matrix))) {
    protein_i <- protein_rows[[i]]
    if (is.na(protein_i)) {
      next
    }
    data <- data.frame(
      id = sample_ids,
      glyco = as.numeric(glyco_matrix[i, sample_ids]),
      protein = as.numeric(protein_matrix[protein_i, sample_ids])
    )
    data <- stats::na.omit(data)
    if (nrow(data) < 3) {
      next
    }
    residuals <- stats::resid(stats::lm(glyco ~ protein, data = data))
    names(residuals) <- data$id
    adjusted[[i]] <- residuals
  }

  keep <- !vapply(adjusted, is.null, logical(1))
  adjusted_matrix <- do.call(rbind, adjusted[keep])
  annotation <- glyco_annotation[match(rownames(adjusted_matrix), glyco_annotation$Gene.Sequence), , drop = FALSE]
  rownames(annotation) <- annotation$Gene.Sequence

  list(matrix = adjusted_matrix, annotation = annotation)
}

#' Extract sex-specific glycopeptide trajectory score matrices
#'
#' This is the package analogue of the manuscript `get_terms(..., 2)` helper
#' used by the glyco ADO survival scripts. It preserves the three returned
#' matrix types: discovery tumor score, reference trajectory score, and
#' discovery trajectory score.
#'
#' @param female_trajectory Female trajectory list, usually `f.gene.mat.adj.list`.
#' @param male_trajectory Male trajectory list, usually `m.gene.mat.adj.list`.
#' @param discovery_clinical Discovery clinical data with `id` and `Gender`.
#' @param reference_clinical Reference clinical data with `id` and `Gender`.
#' @param sex One of `"Male"` or `"Female"`.
#' @param slot Trajectory list element to use. The manuscript glyco SA scripts
#'   use adaptive slot `2`.
#'
#' @return A list with `score.adj`, `dtt.v`, and `dtt.h` matrices.
#' @export
ageTMP_extract_glyco_trajectory_terms <- function(
  female_trajectory,
  male_trajectory,
  discovery_clinical,
  reference_clinical,
  sex = c("Male", "Female"),
  slot = 2
) {
  sex <- match.arg(sex)
  trajectory <- if (identical(sex, "Male")) male_trajectory else female_trajectory
  if (length(trajectory) < slot) {
    stop("Trajectory slot ", slot, " is not available.", call. = FALSE)
  }
  trajectory <- trajectory[[slot]]

  required_matrices <- c("hope.mat.adj", "hope.mat0.adj", "vali.mat0.adj")
  missing <- setdiff(required_matrices, names(trajectory))
  if (length(missing) > 0) {
    stop("Trajectory object is missing matrix element(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  discovery_clinical <- data.frame(discovery_clinical, stringsAsFactors = FALSE)
  reference_clinical <- data.frame(reference_clinical, stringsAsFactors = FALSE)
  discovery_ids <- ageTMP_normalize_sample_ids(discovery_clinical$id)
  reference_ids <- ageTMP_normalize_sample_ids(reference_clinical$id)
  discovery_sex <- standardize_sex(discovery_clinical$Gender)
  reference_sex <- standardize_sex(reference_clinical$Gender)

  select_rows_complete <- function(mat, ids) {
    mat <- as.matrix(mat)
    mat <- mat[, colnames(mat) %in% ids, drop = FALSE]
    t(stats::na.omit(t(mat)))
  }

  list(
    score.adj = select_rows_complete(trajectory$hope.mat.adj, discovery_ids[discovery_sex %in% sex]),
    dtt.v = select_rows_complete(trajectory$vali.mat0.adj, reference_ids[reference_sex %in% sex]),
    dtt.h = select_rows_complete(trajectory$hope.mat0.adj, discovery_ids[discovery_sex %in% sex])
  )
}

#' Prepare discovery glyco ADO-span survival data
#'
#' This wraps the existing discovery survival preparation with manuscript glyco
#' defaults: mutation missingness is zero-filled, age classes are restricted to
#' PED/ADO, and samples are filtered to a supplied glycopeptide sample set.
#'
#' @param clinical Discovery clinical data.
#' @param mutation Optional sample-by-gene mutation indicator matrix.
#' @param glyco_sample_ids Glyco sample IDs available for modeling.
#' @param sex One of `"Male"` or `"Female"`.
#' @param max_day Administrative survival truncation in days.
#'
#' @return A prepared survival data frame.
#' @export
ageTMP_prepare_glyco_ado_discovery_survival <- function(
  clinical,
  mutation = NULL,
  glyco_sample_ids = NULL,
  sex = c("Male", "Female"),
  max_day = 2000
) {
  sex <- match.arg(sex)
  ageTMP_prepare_cdisc_cpsa_survival(
    clinical = clinical,
    mutation = mutation,
    protein_sample_ids = glyco_sample_ids,
    sex = sex,
    age_classes = c("PED", "ADO"),
    max_day = max_day,
    mutation_na = "zero",
    require_complete_mutation = FALSE
  )
}

#' Prepare reference glyco ADO-span survival data
#'
#' @inheritParams ageTMP_prepare_reference_cpsa_survival
#' @param glyco_sample_ids Glyco trajectory sample IDs available for modeling.
#' @param sex One of `"Male"` or `"Female"`.
#'
#' @return A prepared reference survival data frame.
#' @export
ageTMP_prepare_glyco_ado_reference_survival <- function(
  clinical,
  mutation = NULL,
  glyco_sample_ids = NULL,
  sex = c("Male", "Female"),
  max_day = 2000,
  require_complete_mutation = TRUE
) {
  sex <- match.arg(sex)
  ageTMP_prepare_reference_cpsa_survival(
    clinical = clinical,
    mutation = mutation,
    protein_sample_ids = glyco_sample_ids,
    sex = sex,
    age_classes = c("PED", "ADO"),
    max_day = max_day,
    mutation_na = "zero",
    require_complete_mutation = require_complete_mutation
  )
}

#' Fit manuscript glyco ADO-span CPSA models
#'
#' The glyco SA scripts test the combined ADO feature effect using only PED and
#' ADO samples, with PED as the reference age class. This wrapper preserves that
#' model while reusing the package's modality-neutral CPSA engine.
#'
#' @param feature_matrix Numeric glycopeptide-by-sample matrix.
#' @param survival_data Prepared survival data.
#' @param features Optional glycopeptide subset.
#' @param cohort One of `"discovery"` or `"reference"`.
#' @param scale_features Whether to row-scale feature values inside the model.
#'   The manuscript scripts scaled rows immediately before model fitting.
#' @param n_cores Number of cores for fitting.
#' @param progress Whether to print fitting progress.
#' @param progress_every Progress chunk size.
#'
#' @return A CPSA result data frame.
#' @export
ageTMP_fit_glyco_ado_cpsa <- function(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  cohort = c("discovery", "reference"),
  scale_features = TRUE,
  n_cores = 1,
  progress = FALSE,
  progress_every = 500
) {
  cohort <- match.arg(cohort)
  ageTMP_fit_glyco_ado_cpsa_with_spec(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    features = features,
    spec = ageTMP_glyco_ado_cpsa_spec(cohort = cohort, scale_features = scale_features),
    n_cores = n_cores,
    progress = progress,
    progress_every = progress_every
  )
}

#' Manuscript glyco ADO Cox model specification
#'
#' This returns the exact covariate and interaction specification used by the
#' glyco PED/ADO survival scripts. It is useful when callers want the package
#' CPSA engine but need the manuscript glyco ADO settings to be explicit rather
#' than implied by defaults.
#'
#' @param cohort One of `"discovery"` or `"reference"`.
#' @param scale_features Whether features should be row-scaled by the fit
#'   engine. The legacy scripts scaled rows immediately before model fitting.
#'
#' @return A CPSA specification list.
#' @export
ageTMP_glyco_ado_cpsa_spec <- function(
  cohort = c("discovery", "reference"),
  scale_features = TRUE
) {
  cohort <- match.arg(cohort)
  base_covariates <- if (identical(cohort, "reference")) {
    c("age", "mdg", "gbm")
  } else {
    c("age", "Grade", "Cortical", "Midline", "dasc")
  }
  sparse_covariates <- if (identical(cohort, "reference")) {
    c("mdg", "gbm")
  } else {
    c("Grade", "Cortical", "Midline", "dasc")
  }
  sparse_covariates <- c(sparse_covariates, "IDH1_mut", "TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut")

  ageTMP_cpsa_spec(
    feature_col = "score2",
    age_reference = "PED",
    strata = "ADO",
    include_all_combined_test = FALSE,
    base_covariates = base_covariates,
    adjustment_covariates = c("IDH1_mut", "TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut"),
    sparse_covariates = sparse_covariates,
    sparse_min_count = 3,
    sparse_count_method = "factor",
    scale_features = scale_features
  )
}

#' Fit glyco ADO CPSA models with an explicit specification
#'
#' @inheritParams ageTMP_fit_glyco_ado_cpsa
#' @param spec CPSA specification. Defaults to [ageTMP_glyco_ado_cpsa_spec()].
#'
#' @return A CPSA result data frame.
#' @export
ageTMP_fit_glyco_ado_cpsa_with_spec <- function(
  feature_matrix,
  survival_data,
  features = rownames(feature_matrix),
  spec = ageTMP_glyco_ado_cpsa_spec(cohort = "discovery"),
  n_cores = 1,
  progress = FALSE,
  progress_every = 500
) {
  ageTMP_fit_reference_cpsa(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    features = features,
    spec = spec,
    n_cores = n_cores,
    progress = progress,
    progress_every = progress_every
  )
}

#' Add manuscript PED columns to a glyco ADO CPSA result
#'
#' The glyco supplementary table was generated from PED/ADO-only Cox models.
#' The likelihood-ratio p-value is therefore shared by PED and ADO, while the
#' PED sign uses the main feature coefficient and the ADO sign uses the
#' feature-plus-interaction coefficient already present in `ADO.comb.*`.
#'
#' @param fit A glyco CPSA result data frame with `ADO.comb.p`,
#'   `ADO.comb.fdr`, and `score2.coef`.
#'
#' @return `fit` with `PED.comb.*` columns added.
#' @export
ageTMP_add_glyco_ped_columns <- function(fit) {
  fit <- data.frame(fit, check.names = FALSE)
  required <- c("ADO.comb.p", "ADO.comb.fdr", "score2.coef")
  missing <- setdiff(required, names(fit))
  if (length(missing) > 0) {
    stop("Glyco CPSA result is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  fit$PED.comb.p <- fit$ADO.comb.p
  fit$PED.comb.fdr <- fit$ADO.comb.fdr
  fit$PED.comb.coef <- fit$score2.coef
  fit$PED.comb.signed.p <- log10(fit$PED.comb.p) * (-1) * sign(fit$PED.comb.coef)
  fit$PED.comb.signed.fdr <- log10(fit$PED.comb.fdr) * (-1) * sign(fit$PED.comb.coef)
  fit
}

#' Add manuscript glyco local-FDR columns to a validation result
#'
#' This is the package analogue of the glyco `vali.x0` post-processing block in
#' `get_glycosig.R`. It converts a combined p-value column to one-tailed z
#' scores, fits `locfdr`, and adds `.locfdr.sig`, `.locfdr.cont`, and
#' `.signed.locfdr.cont` columns for PED and ADO using the same locfdr threshold.
#'
#' @param fit A glyco validation CPSA result data frame.
#' @param p_column Combined p-value column used for the locfdr fit. The
#'   manuscript no-adjusted glyco blocks use `"ADO.comb.p"` for both PED and ADO.
#' @param classes Combined effect labels to annotate, usually `c("ADO", "PED")`.
#' @param df,nulltype,pct0,bre Arguments passed to [locfdr::locfdr()].
#' @param fdr_cut Local-FDR cutoff used for the discrete `.locfdr.sig` column.
#' @param grid Numeric grid for the continuous local-FDR lookup.
#' @param main Optional plot title passed to `locfdr`.
#'
#' @return `fit` with glyco locfdr columns added.
#' @export
ageTMP_add_glyco_locfdr_columns <- function(
  fit,
  p_column = "ADO.comb.p",
  classes = c("ADO", "PED"),
  df = 40,
  nulltype = 1,
  pct0 = NULL,
  bre = NULL,
  fdr_cut = 0.1,
  grid = seq(from = 0, to = 1, by = 0.0001),
  main = NULL
) {
  if (!requireNamespace("locfdr", quietly = TRUE)) {
    stop("Package `locfdr` is required to add glyco local-FDR columns.", call. = FALSE)
  }
  fit <- ageTMP_add_glyco_ped_columns(fit)
  if (!p_column %in% names(fit)) {
    stop("Glyco CPSA result is missing p-value column: ", p_column, call. = FALSE)
  }

  p <- suppressWarnings(as.numeric(fit[[p_column]]))
  z0 <- stats::qnorm(p, lower.tail = TRUE)
  locfdr_args <- list(zz = z0, df = df, nulltype = nulltype, main = main)
  if (!is.null(pct0)) {
    locfdr_args$pct0 <- pct0
  }
  if (!is.null(bre)) {
    locfdr_args$bre <- bre
  }
  mat <- do.call(locfdr::locfdr, locfdr_args)$mat

  x <- mat[mat[, "Fdrleft"] < fdr_cut, "x"]
  cut <- min(z0, na.rm = TRUE)
  if (length(x) > 0) {
    cut <- max(x, na.rm = TRUE)
  }
  locfdr_sig <- as.numeric(z0 < cut)
  locfdr_cont <- ageTMP_glyco_locfdr_continuous(mat = mat, z0 = z0, grid = grid)

  for (class in classes) {
    signed_col <- paste0(class, ".comb.signed.p")
    if (!signed_col %in% names(fit)) {
      stop("Glyco CPSA result is missing signed p-value column: ", signed_col, call. = FALSE)
    }
    fit[[paste0(class, ".comb.locfdr.sig")]] <- locfdr_sig * sign(fit[[signed_col]])
    fit[[paste0(class, ".comb.locfdr.cont")]] <- locfdr_cont
    fit[[paste0(class, ".comb.signed.locfdr.cont")]] <- log10(locfdr_cont) * (-1) * sign(fit[[signed_col]])
  }

  fit
}

ageTMP_glyco_locfdr_continuous <- function(mat, z0, grid = seq(from = 0, to = 1, by = 0.0001)) {
  locfdr_mat <- vapply(grid, function(y) {
    x <- mat[mat[, "Fdrleft"] < y, "x"]
    cut <- min(z0, na.rm = TRUE)
    if (length(x) > 0) {
      cut <- max(x, na.rm = TRUE)
    }
    as.numeric(z0 < cut)
  }, numeric(length(z0)))

  apply(locfdr_mat, 1, function(x) {
    y <- match(1, x)[1]
    ifelse(is.na(y), 1, grid[[y]])
  })
}

#' Return manuscript glyco local-FDR settings
#'
#' @param sex One of `"Female"` or `"Male"`.
#' @param adjusted Whether the glyco result is protein-adjusted.
#'
#' @return A list of arguments suitable for [ageTMP_add_glyco_locfdr_columns()].
#' @export
ageTMP_glyco_locfdr_settings <- function(
  sex = c("Female", "Male"),
  adjusted = FALSE
) {
  sex <- match.arg(sex)
  if (identical(sex, "Female")) {
    return(list(p_column = "ADO.comb.p", df = 40, nulltype = 1, main = "Female - ADO (glyco)"))
  }
  title <- if (isTRUE(adjusted)) "Male - PED (glyco)" else "Male - ADO (glyco no adj)"
  list(
    p_column = "ADO.comb.p",
    df = 30,
    nulltype = 3,
    pct0_from_z_thresholds = c(-2.2, -1.5),
    main = title
  )
}

#' Add glyco local-FDR columns with manuscript sex-specific settings
#'
#' @inheritParams ageTMP_glyco_locfdr_settings
#' @inheritParams ageTMP_add_glyco_locfdr_columns
#' @param ... Additional arguments passed to [ageTMP_add_glyco_locfdr_columns()].
#'
#' @return `fit` with glyco locfdr columns added.
#' @export
ageTMP_add_glyco_locfdr_by_manuscript_settings <- function(
  fit,
  sex = c("Female", "Male"),
  adjusted = FALSE,
  ...
) {
  settings <- ageTMP_glyco_locfdr_settings(sex = sex, adjusted = adjusted)
  pct0_thresholds <- settings$pct0_from_z_thresholds
  settings$pct0_from_z_thresholds <- NULL
  if (!is.null(pct0_thresholds)) {
    p <- suppressWarnings(as.numeric(ageTMP_add_glyco_ped_columns(fit)[[settings$p_column]]))
    z0 <- stats::qnorm(p, lower.tail = TRUE)
    settings$pct0 <- c(mean(z0 < pct0_thresholds[[1]], na.rm = TRUE), mean(z0 < pct0_thresholds[[2]], na.rm = TRUE))
  }
  do.call(ageTMP_add_glyco_locfdr_columns, c(list(fit = fit), settings, list(...)))
}

#' Assemble manuscript STable4 SA-Glyco-Disc columns from glyco CPSA results
#'
#' @param female_adj,male_adj,female_noadj,male_noadj Glyco CPSA result data
#'   frames, usually the `hope.x00` objects from the four manuscript glyco
#'   ADO result workspaces.
#' @param membrane_features Optional character vector, data frame, or package
#'   membrane-feature table identifying membrane-localized glycopeptides.
#'
#' @return A data frame with the published `SA-Glyco-Disc` columns.
#' @export
ageTMP_build_sa_glyco_disc_from_results <- function(
  female_adj,
  male_adj,
  female_noadj,
  male_noadj,
  membrane_features = NULL
) {
  select_cols <- function(fit) {
    fit <- ageTMP_add_glyco_ped_columns(fit)
    keep <- c(
      "pathway",
      "PED.comb.signed.p",
      "ADO.comb.signed.p",
      "PED.comb.signed.fdr",
      "ADO.comb.signed.fdr"
    )
    missing <- setdiff(keep, names(fit))
    if (length(missing) > 0) {
      stop("Glyco CPSA result is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    fit[, keep, drop = FALSE]
  }

  female_adj <- select_cols(female_adj)
  male_adj <- select_cols(male_adj)
  female_noadj <- select_cols(female_noadj)
  male_noadj <- select_cols(male_noadj)

  adj <- merge(female_adj, male_adj, all = TRUE, by = "pathway", suffixes = c(".female", "male"))
  noadj <- merge(female_noadj, male_noadj, all = TRUE, by = "pathway", suffixes = c(".female", "male"))
  out <- merge(adj, noadj, by = "pathway", all = TRUE, suffixes = c(".adj", ".noadj"))

  if (is.null(membrane_features)) {
    membrane_ids <- character()
  } else if (is.data.frame(membrane_features)) {
    membrane_col <- intersect(c("glycopeptide", "Gene.Sequence", "pathway"), names(membrane_features))
    if (length(membrane_col) == 0) {
      stop("Membrane feature table must contain `glycopeptide`, `Gene.Sequence`, or `pathway`.", call. = FALSE)
    }
    membrane_ids <- membrane_features[[membrane_col[[1]]]]
  } else {
    membrane_ids <- membrane_features
  }

  out$membrane.loc <- ifelse(out$pathway %in% membrane_ids, "Yes", "No")
  names(out)[names(out) == "pathway"] <- "glycoppetide"
  ordered <- c(
    "glycoppetide",
    "membrane.loc",
    "PED.comb.signed.p.female.adj",
    "ADO.comb.signed.p.female.adj",
    "PED.comb.signed.fdr.female.adj",
    "ADO.comb.signed.fdr.female.adj",
    "PED.comb.signed.pmale.adj",
    "ADO.comb.signed.pmale.adj",
    "PED.comb.signed.fdrmale.adj",
    "ADO.comb.signed.fdrmale.adj",
    "PED.comb.signed.p.female.noadj",
    "ADO.comb.signed.p.female.noadj",
    "PED.comb.signed.fdr.female.noadj",
    "ADO.comb.signed.fdr.female.noadj",
    "PED.comb.signed.pmale.noadj",
    "ADO.comb.signed.pmale.noadj",
    "PED.comb.signed.fdrmale.noadj",
    "ADO.comb.signed.fdrmale.noadj"
  )
  out[, ordered, drop = FALSE]
}
