\name{ageTMP_fit_reference_protein_cpsa}
\alias{ageTMP_fit_reference_protein_cpsa}
\title{Fit reference-cohort protein CPSA models}
\usage{
ageTMP_fit_reference_protein_cpsa(
  protein_matrix,
  survival_data,
  genes = rownames(protein_matrix),
  age_reference = "ADULT",
  mutation_genes = c("IDH1", "TP53", "ATRX", "H33A", "ATM")
)
}
\arguments{
\item{protein_matrix}{Numeric gene-by-sample protein matrix.}

\item{survival_data}{Prepared survival covariates from
\code{\link[=ageTMP_prepare_reference_cpsa_survival]{ageTMP_prepare_reference_cpsa_survival()}}.}

\item{genes}{Optional subset of genes to fit.}

\item{age_reference}{Reference age class for the interaction model.}

\item{mutation_genes}{Mutation gene names without \code{_mut}.}
}
\value{
A data frame with coefficients, p-values, BY-adjusted FDR values, and signed
combined age-stratum statistics.
}
\description{
Protein convenience alias for the modality-neutral
\code{\link[=ageTMP_fit_reference_cpsa]{ageTMP_fit_reference_cpsa()}} engine.
}
