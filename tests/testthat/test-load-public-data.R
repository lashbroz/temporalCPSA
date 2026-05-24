test_that("data source manifest is well formed", {
  sources <- ageTMP_data_sources(data_dir = tempdir())
  expect_true(all(c("source", "file", "path", "exists") %in% names(sources)))
  expect_true("STable4" %in% sources$source)
  expect_false(any(sources$exists))
})

test_that("sample ID normalization matches paper helper behavior", {
  expect_equal(
    ageTMP_normalize_sample_ids(c("X7316.1052", "A.7316.1099", "P.7316.114")),
    c("7316-1052", "7316-1099", "7316-114")
  )
})

test_that("feature matrix loader supports multiple modalities", {
  data_dir <- Sys.getenv("AGETMP_TEST_DATA_DIR", unset = "data")
  skip_if_not(file.exists(file.path(data_dir, "cDisc_proteome_imputed_data_09152023.tsv")))
  skip_if_not(file.exists(file.path(data_dir, "cDisc_rna_coding_10192023.tsv")))

  matrices <- ageTMP_load_feature_matrix(
    data_dir = data_dir,
    modality = c("protein", "rna")
  )

  expect_named(matrices, c("protein", "rna"))
  expect_true(all(vapply(matrices, is.matrix, logical(1))))
  expect_true(all(vapply(matrices, nrow, integer(1)) > 0))
})

test_that("phospho feature matrix defaults are explicit by aggregation level", {
  data_dir <- Sys.getenv("AGETMP_TEST_DATA_DIR", unset = "data")
  skip_if_not(file.exists(file.path(data_dir, "cDisc_phosphosite_imputed_data_ischemia_removed_motif_11032023.tsv")))

  site_matrix <- ageTMP_load_feature_matrix(
    data_dir = data_dir,
    modality = "phospho",
    collapse = FALSE
  )
  gene_matrix <- ageTMP_load_feature_matrix(
    data_dir = data_dir,
    modality = "phospho",
    collapse = TRUE
  )

  expect_true(any(grepl("^NP_", rownames(site_matrix))))
  expect_identical(anyDuplicated(rownames(gene_matrix)), 0L)
  expect_true(nrow(site_matrix) > nrow(gene_matrix))
})
