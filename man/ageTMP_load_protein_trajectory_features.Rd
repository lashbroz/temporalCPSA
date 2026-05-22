\name{ageTMP_load_protein_trajectory_features}
\alias{ageTMP_load_protein_trajectory_features}
\title{Load manuscript protein trajectory feature universe}
\usage{
ageTMP_load_protein_trajectory_features(
  rdata_path,
  list_index = 2,
  matrix_name = "hope.mat0.adj"
)
}
\arguments{
\item{rdata_path}{Path to a \code{protein_rev_tadj62_list.RData}-style object.}

\item{list_index}{List element containing the manuscript trajectory matrices.}

\item{matrix_name}{Matrix whose row names define the protein feature universe.}
}
\value{
Character vector of protein features.
}
\description{
Load manuscript protein trajectory feature universe
}
