test_that("trajectory SD ranking and filtering are explicit", {
  trajectory <- data.frame(
    feature = rep(c("flat", "dynamic"), each = 6),
    sex = rep(rep(c("Male", "Female"), each = 3), times = 2),
    tissue = "Tumor",
    fit = c(
      0.01, 0.02, 0.03,
      0.01, 0.02, 0.03,
      -1, 0, 1,
      -0.8, 0, 0.8
    )
  )

  sd_rank <- ageTMP_rank_trajectory_sd(
    trajectory,
    group_cols = "sex",
    tissue = "Tumor"
  )

  expect_true(all(c("feature", "sex", "sd", "rank_sd") %in% names(sd_rank)))
  expect_true(all(sd_rank$sd[sd_rank$feature == "dynamic"] > 0.15))
  expect_true(all(sd_rank$sd[sd_rank$feature == "flat"] < 0.15))

  kept <- ageTMP_filter_trajectory_sd(sd_rank, sd_min = 0.15, group_cols = "sex")
  expect_equal(kept, "dynamic")

  kept_per_group <- ageTMP_filter_trajectory_sd(
    sd_rank,
    sd_min = 0.15,
    group_cols = "sex",
    keep = "per_group"
  )
  expect_equal(unique(kept_per_group$feature), "dynamic")
})
