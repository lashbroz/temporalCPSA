test_that("age class derivation follows manuscript five-class intervals", {
  age <- c(0, 15, 15.1, 26, 30, 40, 62, 63, 80, 81, NA)
  classes <- ageTMP_derive_age_class(age)
  ranges <- ageTMP_derive_age_class_range(age)

  expect_equal(
    as.character(classes),
    c("PED", "PED", "ADO", "ADO", "YA", "YA", "ADULT", "SEN", "SEN", NA, NA)
  )
  expect_equal(
    as.character(ranges),
    c("[0,15]", "[0,15]", "(15,26]", "(15,26]", "(26,40]", "(26,40]", "(40,62]", "(62,80]", "(62,80]", NA, NA)
  )
})

test_that("age class columns can be added to a clinical data frame", {
  clinical <- data.frame(
    id = c("S1", "S2", "S3"),
    age = c(5, 20, 70),
    stringsAsFactors = FALSE
  )

  out <- ageTMP_add_age_class(clinical)

  expect_true(all(c("age_class", "age_class_range") %in% names(out)))
  expect_equal(as.character(out$age_class), c("PED", "ADO", "SEN"))
  expect_equal(as.character(out$age_class_range), c("[0,15]", "(15,26]", "(62,80]"))
})
