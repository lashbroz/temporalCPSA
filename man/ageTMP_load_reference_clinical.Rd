\name{ageTMP_load_reference_clinical}
\alias{ageTMP_load_reference_clinical}
\title{Load reference-cohort clinical data}
\usage{
ageTMP_load_reference_clinical(data_dir = "data")
}
\arguments{
\item{data_dir}{Path to the public data directory.}
}
\value{
A clinical data frame.
}
\description{
Public-facing alias for \code{\link[=ageTMP_load_cdisc_clinical]{ageTMP_load_cdisc_clinical()}}. In the manuscript
source files, the external/reference survival cohort is labeled \code{cDisc}; in
the package API, "reference cohort" is the preferred general term.
}
