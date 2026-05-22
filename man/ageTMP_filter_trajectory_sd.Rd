\name{ageTMP_filter_trajectory_sd}
\alias{ageTMP_filter_trajectory_sd}
\title{Filter trajectory features by standard deviation}
\usage{
ageTMP_filter_trajectory_sd(
  sd_rank,
  sd_min = 0.15,
  group_cols = c("sex"),
  keep = c("all_groups", "any_group", "per_group")
)
}
\arguments{
\item{sd_rank}{Data frame returned by \code{\link{ageTMP_rank_trajectory_sd}}.}
\item{sd_min}{Minimum standard deviation required to keep a feature.}
\item{group_cols}{Columns defining groups that must pass the threshold.}
\item{keep}{One of \code{"all_groups"}, \code{"any_group"}, or \code{"per_group"}.}
}
\value{
A character vector of kept features for \code{"all_groups"} or \code{"any_group"},
and a filtered data frame for \code{"per_group"}.
}
\description{
Apply an explicit dynamic-range filter to trajectory SD rankings. The default
keeps features with \code{sd >= 0.15}, matching the omit-list threshold used in
the manuscript \code{rank_df_rev.RData} workflow.
}
