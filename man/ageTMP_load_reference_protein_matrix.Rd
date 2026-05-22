\name{ageTMP_load_reference_protein_matrix}
\alias{ageTMP_load_reference_protein_matrix}
\title{Load reference-cohort protein matrix}
\usage{
ageTMP_load_reference_protein_matrix(data_dir = "data", collapse = TRUE)
}
\arguments{
\item{data_dir}{Path to the public data directory.}

\item{collapse}{Whether to average rows with the same gene symbol.}
}
\value{
A numeric gene-by-sample protein matrix.
}
\description{
Public-facing alias for \code{\link[=ageTMP_load_cdisc_protein_matrix]{ageTMP_load_cdisc_protein_matrix()}}. The manuscript
implementation currently supports protein CPSA first; the naming leaves room
for RNA, glyco, phosphosite, and other feature matrices to use the same
reference-cohort vocabulary.
}
