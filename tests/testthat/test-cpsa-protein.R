test_that("cDisc mutation loader returns requested mutation indicators", {
  data_dir <- file.path(getwd(), "..", "data")
  skip_if_not(file.exists(file.path(data_dir, "cDisc_mutation_10192023.tsv")))

  mutation <- ageTMP_load_cdisc_mutation(data_dir = data_dir)

  expect_true(all(c("IDH1_mut", "TP53_mut", "ATRX_mut", "H33A_mut", "ATM_mut") %in% colnames(mutation)))
  expect_true(is.numeric(mutation[, "TP53_mut"]))
  expect_true(nrow(mutation) >= 188)
})

test_that("cDisc clinical survival preparation applies manuscript filters", {
  data_dir <- file.path(getwd(), "..", "data")
  skip_if_not(file.exists(file.path(data_dir, "STable1.xlsx")))

  clinical <- ageTMP_load_cdisc_clinical(data_dir = data_dir)
  mutation <- ageTMP_load_cdisc_mutation(data_dir = data_dir)
  protein <- ageTMP_load_cdisc_protein_matrix(data_dir = data_dir)
  surv <- ageTMP_prepare_cdisc_cpsa_survival(
    clinical = clinical,
    mutation = mutation,
    protein_sample_ids = colnames(protein),
    sex = "Male"
  )

  expect_true(all(surv$Gender == "Male"))
  expect_true(all(surv$treat_status == "Treatment naive"))
  expect_true(all(surv$First.Diagnosis == "Yes"))
  expect_true(all(!is.na(surv$days)))
})

test_that("cDisc survival preparation can preserve mutation missingness", {
  clinical <- data.frame(
    id = c("S1", "S2"),
    cDisc_os = c(10, 20),
    cDisc_os_status = c(1, 0),
    cDisc_pfs = c(10, 20),
    cDisc_pfs_status = c(1, 0),
    cDisc_clinical_status_at_collection_event = c("Alive", "Alive"),
    cDisc_treat_status = c("Treatment naive", "Treatment naive"),
    cDisc_First_Diagnosis = c("Yes", "Yes"),
    cDisc_age = c(12, 16),
    cDisc_Gender = c("Female", "Female"),
    cDisc_WHO_Grade = c(3, 3),
    cDisc_tumor_loc = c("Cortical", "Midline"),
    stringsAsFactors = FALSE
  )
  mutation <- matrix(NA_real_, nrow = 1, ncol = 1)
  rownames(mutation) <- "S1"
  colnames(mutation) <- "TP53_mut"

  preserved <- ageTMP_prepare_cdisc_cpsa_survival(
    clinical = clinical,
    mutation = mutation,
    mutation_na = "preserve"
  )
  zeroed <- ageTMP_prepare_cdisc_cpsa_survival(
    clinical = clinical,
    mutation = mutation,
    mutation_na = "zero"
  )

  expect_true(any(is.na(preserved$TP53_mut)))
  expect_false(any(is.na(zeroed$TP53_mut)))
})

test_that("generic CPSA can fit without age-class interactions", {
  set.seed(1)
  sample_ids <- paste0("S", seq_len(40))
  feature_matrix <- matrix(
    rnorm(80),
    nrow = 2,
    dimnames = list(c("feature_a", "feature_b"), sample_ids)
  )
  survival_data <- data.frame(
    days = rexp(40, 0.02) + 1,
    os.status = c(rep(1, 5), rep(0, 5), rbinom(30, 1, 0.6)),
    age = rnorm(40, 30, 10),
    row.names = sample_ids
  )

  fit <- ageTMP_fit_reference_cpsa(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    spec = ageTMP_cpsa_spec(
      age_class_col = NULL,
      age_reference = NULL,
      strata = NULL,
      include_all_combined_test = FALSE,
      base_covariates = "age",
      adjustment_covariates = character(),
      remove_sparse_covariates = FALSE,
      scale_features = FALSE
    )
  )

  expect_equal(nrow(fit), 2)
  expect_true(all(c("pathway", "feature_score.coef", "feature_score.p", "feature_score.fdr") %in% names(fit)))
  expect_false(any(grepl("[.]comb[.]", names(fit))))
})

test_that("generic CPSA can screen features by trajectory SD", {
  set.seed(11)
  sample_ids <- paste0("S", seq_len(40))
  feature_matrix <- matrix(
    rnorm(120),
    nrow = 3,
    dimnames = list(c("flat", "dynamic", "other"), sample_ids)
  )
  survival_data <- data.frame(
    days = rexp(40, 0.02) + 1,
    os.status = c(rep(1, 5), rep(0, 5), rbinom(30, 1, 0.6)),
    age = rnorm(40, 30, 10),
    row.names = sample_ids
  )
  spec <- ageTMP_cpsa_spec(
    age_class_col = NULL,
    age_reference = NULL,
    strata = NULL,
    include_all_combined_test = FALSE,
    base_covariates = "age",
    adjustment_covariates = character(),
    remove_sparse_covariates = FALSE,
    scale_features = FALSE
  )

  fit <- ageTMP_fit_reference_cpsa(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    spec = spec,
    trajectory_sd = c(flat = 0.03, dynamic = 0.22, other = 0.12),
    trajectory_sd_min = 0.15
  )

  expect_equal(fit$pathway, "dynamic")
})

test_that("generic CPSA can screen grouped trajectory SD data frames", {
  set.seed(12)
  sample_ids <- paste0("S", seq_len(40))
  feature_matrix <- matrix(
    rnorm(80),
    nrow = 2,
    dimnames = list(c("passes_one", "passes_both"), sample_ids)
  )
  survival_data <- data.frame(
    days = rexp(40, 0.02) + 1,
    os.status = c(rep(1, 5), rep(0, 5), rbinom(30, 1, 0.6)),
    age = rnorm(40, 30, 10),
    row.names = sample_ids
  )
  sd_rank <- data.frame(
    sex = rep(c("Male", "Female"), each = 2),
    feature = rep(c("passes_one", "passes_both"), times = 2),
    sd = c(0.20, 0.21, 0.05, 0.22)
  )
  spec <- ageTMP_cpsa_spec(
    age_class_col = NULL,
    age_reference = NULL,
    strata = NULL,
    include_all_combined_test = FALSE,
    base_covariates = "age",
    adjustment_covariates = character(),
    remove_sparse_covariates = FALSE,
    scale_features = FALSE
  )

  any_group <- ageTMP_fit_reference_cpsa(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    spec = spec,
    trajectory_sd = sd_rank,
    trajectory_sd_min = 0.15,
    trajectory_sd_group_cols = "sex",
    trajectory_sd_keep = "any_group"
  )
  all_groups <- ageTMP_fit_reference_cpsa(
    feature_matrix = feature_matrix,
    survival_data = survival_data,
    spec = spec,
    trajectory_sd = sd_rank,
    trajectory_sd_min = 0.15,
    trajectory_sd_group_cols = "sex",
    trajectory_sd_keep = "all_groups"
  )

  expect_equal(sort(any_group$pathway), c("passes_both", "passes_one"))
  expect_equal(all_groups$pathway, "passes_both")
})

test_that("matrix Cox engine matches formula Cox engine without age classes", {
  set.seed(2)
  sample_ids <- paste0("S", seq_len(80))
  feature_matrix <- matrix(
    rnorm(160),
    nrow = 2,
    dimnames = list(c("feature_a", "feature_b"), sample_ids)
  )
  survival_data <- data.frame(
    days = rexp(80, 0.02) + 1,
    os.status = rbinom(80, 1, 0.65),
    age = rnorm(80, 30, 10),
    row.names = sample_ids
  )
  base_spec <- ageTMP_cpsa_spec(
    age_class_col = NULL,
    age_reference = NULL,
    strata = NULL,
    include_all_combined_test = FALSE,
    base_covariates = "age",
    adjustment_covariates = character(),
    remove_sparse_covariates = FALSE,
    scale_features = FALSE
  )

  formula_fit <- ageTMP_fit_reference_cpsa(feature_matrix, survival_data, spec = base_spec)
  matrix_fit <- ageTMP_fit_reference_cpsa(
    feature_matrix,
    survival_data,
    spec = utils::modifyList(base_spec, list(engine = "coxph.fit"))
  )

  numeric_cols <- names(formula_fit)[vapply(formula_fit, is.numeric, logical(1))]
  expect_equal(matrix_fit[, numeric_cols], formula_fit[, numeric_cols], tolerance = 1e-8)
})

test_that("matrix Cox engine matches formula Cox engine with age-class interactions", {
  set.seed(3)
  sample_ids <- paste0("S", seq_len(120))
  feature_matrix <- matrix(
    rnorm(240),
    nrow = 2,
    dimnames = list(c("feature_a", "feature_b"), sample_ids)
  )
  survival_data <- data.frame(
    days = rexp(120, 0.018) + 1,
    os.status = rbinom(120, 1, 0.7),
    age = rnorm(120, 28, 9),
    age_class_new = rep(c("PED", "ADO", "YA", "ADULT"), each = 30),
    row.names = sample_ids
  )
  base_spec <- ageTMP_cpsa_spec(
    base_covariates = "age",
    adjustment_covariates = character(),
    remove_sparse_covariates = FALSE,
    scale_features = FALSE
  )

  formula_fit <- ageTMP_fit_reference_cpsa(feature_matrix, survival_data, spec = base_spec)
  matrix_fit <- ageTMP_fit_reference_cpsa(
    feature_matrix,
    survival_data,
    spec = utils::modifyList(base_spec, list(engine = "coxph.fit"))
  )

  numeric_cols <- names(formula_fit)[vapply(formula_fit, is.numeric, logical(1))]
  expect_equal(matrix_fit[, numeric_cols], formula_fit[, numeric_cols], tolerance = 1e-8)
})
