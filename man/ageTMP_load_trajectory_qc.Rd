\name{ageTMP_load_trajectory_qc}
\alias{ageTMP_load_trajectory_qc}
\title{Load manuscript trajectory QC omit lists}
\usage{
ageTMP_load_trajectory_qc(rdata_path, type = "protein")
}
\arguments{
\item{rdata_path}{Path to a \code{rank_df_rev.RData}-style object.}

\item{type}{Molecular data type to keep, such as \code{"protein"}.}
}
\value{
A list with the rank data frame and male/female omit vectors.
}
\description{
The manuscript CPSA workflow uses sex-specific omit lists for features with
very flat age-dependent trajectory signal. These omit lists are derived from
fitted trajectory matrices by ranking the row-wise standard deviation and
using \code{sd < 0.15} as the default low-dynamic-range threshold.
}
