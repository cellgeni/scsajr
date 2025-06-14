% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scsajr.R
\name{calc_psi}
\alias{calc_psi}
\title{Calculate percent spliced‐in (PSI) per segment}
\usage{
calc_psi(se, min_cov = 10)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} with assays named \code{i} (inclusion counts) and \code{e} (exclusion counts).
Rows are segments; columns are pseudobulk samples.}

\item{min_cov}{Integer; minimum total junction coverage (\code{i + e}) required to compute a valid PSI.
For any \code{i + e < min_cov}, the PSI is set to \code{NA}. Default is \code{10}.}
}
\value{
A numeric matrix of dimensions \code{nrow(se)} × \code{ncol(se)} where each entry is the PSI value
for that segment and sample, or \code{NA} if coverage is insufficient.
}
\description{
Given a \code{SummarizedExperiment} containing inclusion (\code{i}) and exclusion (\code{e}) counts,
this function computes, for each segment and sample, the percent spliced‐in: \code{PSI = i / (i + e)}
If the total counts (\code{i + e}) for a segment in a sample are below \code{min_cov}, that PSI is set to \code{NA}.
}
\examples{
\dontrun{
# Assume 'pseudobulk_se' is a SummarizedExperiment with assays 'i' and 'e':
psi_mat <- calc_psi(pseudobulk_se, min_cov = 10)
head(psi_mat)
}

}
\seealso{
\code{\link{pseudobulk}}, \code{\link{calc_cpm}}
}
