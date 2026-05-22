\name{ageTMP_build_sa_protein_reference}
\alias{ageTMP_build_sa_protein_reference}
\title{Build manuscript STable4 protein reference-cohort CPSA columns}
\usage{
ageTMP_build_sa_protein_reference(
  data_dir = "data",
  genes = NULL,
  mode = c("manuscript", "standardized")
)
}
\arguments{
\item{data_dir}{Path to the public data directory.}

\item{genes}{Optional gene subset for testing.}

\item{mode}{Reproducibility mode. \code{"manuscript"} preserves the original
sex-specific mutation missingness behavior used by the manuscript protein CPSA
scripts: male mutation NAs are set to zero and female mutation NAs are
preserved. \code{"standardized"} zero-fills missing mutation calls for both
sexes.}
}
\value{
A data frame with manuscript reference-cohort signed log10 p-value columns for
male and female.
}
\description{
Public-facing alias for \code{\link[=ageTMP_build_sa_protein_cdisc]{ageTMP_build_sa_protein_cdisc()}}. It returns the
manuscript \code{SA-Protein-cDisc-Ref} reference-cohort signed log10 p-value
columns, preserving the published column naming convention.
}
