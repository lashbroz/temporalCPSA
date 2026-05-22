test_that("normal reference loader fails clearly when package data are absent", {
  expect_error(
    ageTMP_load_normal_reference(path = file.path(tempdir(), "missing.rds")),
    "Normal-reference data were not found"
  )
})

test_that("normal reference includes manuscript protein N-TMP clusters", {
  normal_reference <- ageTMP_load_normal_reference()

  expect_named(normal_reference$protein, c("matrix", "sample_metadata", "clusters"))
  expect_named(
    normal_reference$protein$clusters,
    c("consensus_k4_raw", "consensus_k4", "remap", "source_file", "source_object")
  )
  expect_true(all(c("Male,ALDH9A1", "Female,L1CAM") %in% names(normal_reference$protein$clusters$consensus_k4)))
  expect_equal(as.character(normal_reference$protein$clusters$consensus_k4["Male,ALDH9A1"]), "2")
  expect_equal(as.character(normal_reference$protein$clusters$consensus_k4["Female,L1CAM"]), "1")
})
