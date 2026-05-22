test_that("glyco loader preserves glycopeptide-level source rows", {
  data_dir <- Sys.getenv("AGETMP_TEST_DATA_DIR", unset = "data")
  if (!file.exists(file.path(data_dir, "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv"))) {
    skip("public manuscript data directory is not available")
  }

  glyco <- ageTMP_load_cdisc_glyco_matrix(data_dir = data_dir)

  expect_equal(nrow(glyco$matrix), 3440)
  expect_equal(ncol(glyco$matrix), 77)
  expect_true(all(c("Modification", "Gene", "Gene.Sequence", "GlycanType") %in% names(glyco$annotation)))
  expect_equal(rownames(glyco$matrix), glyco$annotation$Gene.Sequence)
})

test_that("glyco ADO CPSA wrapper can fit a small discovery model", {
  data_dir <- Sys.getenv("AGETMP_TEST_DATA_DIR", unset = "data")
  if (!file.exists(file.path(data_dir, "Disc_glyco_v2_imputed_batch1+2_05082024_011524.tsv"))) {
    skip("public manuscript data directory is not available")
  }

  glyco <- ageTMP_load_cdisc_glyco_matrix(data_dir = data_dir)
  clinical <- ageTMP_load_cdisc_clinical(data_dir = data_dir)
  mutation <- ageTMP_load_cdisc_mutation(data_dir = data_dir)
  surv <- ageTMP_prepare_glyco_ado_discovery_survival(
    clinical = clinical,
    mutation = mutation,
    glyco_sample_ids = colnames(glyco$matrix),
    sex = "Male"
  )

  fit <- ageTMP_fit_glyco_ado_cpsa(
    feature_matrix = glyco$matrix,
    survival_data = surv,
    features = rownames(glyco$matrix)[1:2]
  )

  expect_equal(nrow(fit), 2)
  expect_true(all(c("pathway", "ADO.comb.p", "ADO.comb.signed.p") %in% names(fit)))
})
