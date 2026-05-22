#' Derive age-class labels from age
#'
#' Convert numeric ages to ordered age-class labels using explicit breakpoints.
#' The default breakpoints reproduce the manuscript five-class age grouping:
#' PED, ADO, YA, ADULT, and SEN.
#'
#' @param age Numeric age vector.
#' @param breaks Numeric breakpoints passed to [base::cut()].
#' @param labels Age-class labels. Length must be `length(breaks) - 1`.
#' @param right Whether intervals are closed on the right.
#' @param include_lowest Whether the lowest age should be included in the first
#'   interval.
#' @param ordered Whether to return an ordered factor.
#'
#' @return An ordered factor of age-class labels.
#' @export
ageTMP_derive_age_class <- function(
  age,
  breaks = c(0, 15, 26, 40, 62, 80),
  labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
) {
  if (length(labels) != length(breaks) - 1) {
    stop("`labels` must have length `length(breaks) - 1`.", call. = FALSE)
  }
  age <- suppressWarnings(as.numeric(age))
  cut(
    age,
    breaks = breaks,
    labels = labels,
    right = right,
    include.lowest = include_lowest,
    ordered_result = ordered
  )
}

#' Derive age-class interval labels from age
#'
#' Return the interval labels corresponding to [ageTMP_derive_age_class()]. This
#' is useful when reproducing manuscript tables that show both a semantic class
#' label, such as `PED`, and the numeric range, such as `[0,15]`.
#'
#' @inheritParams ageTMP_derive_age_class
#'
#' @return A factor of interval labels.
#' @export
ageTMP_derive_age_class_range <- function(
  age,
  breaks = c(0, 15, 26, 40, 62, 80),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
) {
  age <- suppressWarnings(as.numeric(age))
  cut(
    age,
    breaks = breaks,
    right = right,
    include.lowest = include_lowest,
    ordered_result = ordered
  )
}

#' Add age-class columns to a clinical data frame
#'
#' @param clinical Clinical data frame.
#' @param age_col Column containing numeric age.
#' @param class_col Name of the derived age-class column to create.
#' @param range_col Optional name of the derived age-range column to create. Set
#'   to `NULL` to omit the range column.
#' @inheritParams ageTMP_derive_age_class
#'
#' @return `clinical` with derived age-class columns appended.
#' @export
ageTMP_add_age_class <- function(
  clinical,
  age_col = "age",
  class_col = "age_class",
  range_col = "age_class_range",
  breaks = c(0, 15, 26, 40, 62, 80),
  labels = c("PED", "ADO", "YA", "ADULT", "SEN"),
  right = TRUE,
  include_lowest = TRUE,
  ordered = TRUE
) {
  if (!age_col %in% names(clinical)) {
    stop("Age column not found: ", age_col, call. = FALSE)
  }
  clinical[[class_col]] <- ageTMP_derive_age_class(
    clinical[[age_col]],
    breaks = breaks,
    labels = labels,
    right = right,
    include_lowest = include_lowest,
    ordered = ordered
  )
  if (!is.null(range_col)) {
    clinical[[range_col]] <- ageTMP_derive_age_class_range(
      clinical[[age_col]],
      breaks = breaks,
      right = right,
      include_lowest = include_lowest,
      ordered = ordered
    )
  }
  clinical
}
