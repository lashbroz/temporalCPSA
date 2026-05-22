test_that("age class supplement builder derives labels and ranges", {
  clinical <- data.frame(
    id = paste0("S", 1:5),
    age = c(1, 16, 30, 50, 70)
  )

  out <- ageTMP_build_age_class_supplement(clinical)

  expect_equal(names(out), c("id", "age.class", "age.class.range", "sil_width"))
  expect_equal(out$age.class, c("PED", "ADO", "YA", "ADULT", "SEN"))
  expect_true(all(is.na(out$sil_width)))
})

test_that("age class supplement builder computes silhouette widths from TMP matrix", {
  clinical <- data.frame(
    id = paste0("S", 1:6),
    age = c(1, 2, 16, 17, 45, 46)
  )
  tmp <- matrix(
    c(
      -2, -1.8, 0, 0.2, 2, 2.2,
      -1, -1.2, 0.1, 0, 1.1, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("feature1", "feature2"), clinical$id)
  )

  out <- ageTMP_build_age_class_supplement(clinical, tmp_matrix = tmp)

  expect_false(any(is.na(out$sil_width)))
  expect_equal(nrow(out), 6)
})
