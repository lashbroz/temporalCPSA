\name{ageTMP_load_reference_mutation}
\alias{ageTMP_load_reference_mutation}
\title{Load reference-cohort mutation covariates}
\usage{
ageTMP_load_reference_mutation(
  data_dir = "data",
  genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
)
}
\arguments{
\item{data_dir}{Path to the public data directory.}

\item{genes}{Character vector of genes to return.}
}
\value{
A numeric sample-by-gene mutation indicator matrix.
}
\description{
Public-facing alias for \code{\link[=ageTMP_load_cdisc_mutation]{ageTMP_load_cdisc_mutation()}}. The returned matrix is
suitable for reference-cohort CPSA adjustment covariates.
}
