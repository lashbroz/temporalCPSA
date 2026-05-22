\name{ageTMP_rank_trajectory_sd}
\alias{ageTMP_rank_trajectory_sd}
\title{Rank trajectory features by standard deviation}
\usage{
ageTMP_rank_trajectory_sd(
  trajectory,
  feature_col = "feature",
  value_col = "fit",
  group_cols = c("sex"),
  tissue_col = "tissue",
  tissue = NULL
)
}
\arguments{
\item{trajectory}{Long-format trajectory data frame.}
\item{feature_col}{Column containing feature names.}
\item{value_col}{Column containing trajectory values, usually \code{fit}.}
\item{group_cols}{Columns defining independent rankings.}
\item{tissue_col}{Optional tissue column used with \code{tissue}.}
\item{tissue}{Optional tissue value to retain before ranking, such as \code{"Tumor"}.}
}
\value{
A data frame with grouping columns, \code{feature}, \code{sd}, and \code{rank_sd}.
}
\description{
Compute within-group standard deviation of fitted trajectory values. This makes
explicit the manuscript QC idea behind \code{rank_df_rev.RData}, where some
downstream analyses exclude features with very flat age-dependent tumor
trajectories, commonly using \code{sd < 0.15} as the omit threshold.
}
