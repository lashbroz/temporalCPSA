test_that("TMP matrices can be selected, combined, and clustered", {
  set.seed(21)
  samples <- paste0("S", seq_len(12))
  mat1 <- matrix(rnorm(60), nrow = 5, dimnames = list(paste0("g", 1:5), samples))
  mat2 <- matrix(rnorm(72), nrow = 6, dimnames = list(paste0("g", 1:6), samples))
  mat1[1, ] <- seq(-2, 2, length.out = length(samples))
  mat2[2, ] <- rep(c(-1, 1), each = 6)

  combined <- ageTMP_combine_tmp_matrices(
    list(rna = mat1, phospho = mat2),
    proportion = 0.5,
    max_features = 3
  )

  expect_equal(ncol(combined), length(samples))
  expect_true(all(grepl("::", rownames(combined))))
  expect_lte(nrow(combined), 6)

  clustered <- ageTMP_cluster_tmp_samples(
    combined,
    k_max = 3,
    reps = 5,
    plot = NULL,
    use_consensus = FALSE
  )

  expect_true(all(c("k2", "k3") %in% names(clustered$clusters)))
  expect_equal(names(clustered$clusters$k2), samples)
})

test_that("TMP age clusters can be summarized by age range", {
  clusters <- c(S1 = 1, S2 = 1, S3 = 2, S4 = 2)
  clinical <- data.frame(
    id = c("S3", "S1", "S4", "S2"),
    age = c(20, 5, 30, 10)
  )

  summary <- ageTMP_summarize_tmp_age_clusters(clusters, clinical)

  expect_equal(as.character(summary$clusters), c("1", "1", "2", "2"))
  expect_equal(nrow(summary$age_range), 2)
})

test_that("AD-TMP age segmentation returns contiguous age intervals", {
  set.seed(42)
  samples <- paste0("S", seq_len(30))
  age <- seq(1, 75, length.out = length(samples))
  mat <- rbind(
    early = c(rep(2, 10), rep(0, 20)),
    mid = c(rep(0, 10), rep(2, 10), rep(0, 10)),
    late = c(rep(0, 20), rep(2, 10))
  )
  mat <- mat + matrix(rnorm(length(mat), sd = 0.02), nrow = nrow(mat))
  dimnames(mat) <- list(rownames(mat), samples)
  names(age) <- samples

  fit <- ageTMP_segment_age_classes(mat, age, k = 3, min_n = 5)

  expect_s3_class(fit, "ageTMP_age_segmentation")
  expect_equal(nrow(fit$segment_summary), 3)
  expect_true(all(diff(fit$ordered_samples$age) >= 0))
  expect_true(all(fit$segment_summary$start_index[-1] == head(fit$segment_summary$end_index, -1) + 1))
  expect_equal(names(fit$classes), fit$ordered_samples$id)
})

test_that("AD-TMP segmentation diagnostic matrix preserves ordered samples", {
  set.seed(84)
  samples <- paste0("S", seq_len(24))
  age <- seq(2, 70, length.out = length(samples))
  mat <- matrix(rnorm(6 * 24), nrow = 6, dimnames = list(paste0("f", 1:6), samples))
  mat[, 1:8] <- mat[, 1:8] - 1
  mat[, 17:24] <- mat[, 17:24] + 1
  names(age) <- samples

  diagnostic <- ageTMP_segment_diagnostic_matrix(mat, age, k_values = 2:4, min_n = 4)

  expect_equal(dim(diagnostic$matrix), c(3, 24))
  expect_equal(rownames(diagnostic$matrix), paste0("k", 2:4))
  expect_equal(colnames(diagnostic$matrix), diagnostic$ordered_samples$id)
  expect_true(all(diff(diagnostic$ordered_samples$age) >= 0))
})

test_that("AD-TMP segmentation diagnostic can overlay suggested age classes", {
  set.seed(126)
  samples <- paste0("S", seq_len(24))
  age <- seq(2, 70, length.out = length(samples))
  mat <- matrix(rnorm(6 * 24), nrow = 6, dimnames = list(paste0("f", 1:6), samples))
  names(age) <- samples
  diagnostic <- ageTMP_segment_diagnostic_matrix(mat, age, k_values = 2:4, min_n = 4)

  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = 600, height = 360)
  expect_invisible(ageTMP_plot_segment_diagnostic(
    diagnostic,
    selected_k = 3,
    suggested_cutpoints = c(15, 40),
    suggested_labels = c("PED", "MID", "ADULT")
  ))
  grDevices::dev.off()

  expect_true(file.exists(file))
})
