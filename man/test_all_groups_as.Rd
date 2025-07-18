% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scsajr.R
\name{test_all_groups_as}
\alias{test_all_groups_as}
\title{Test alternative splicing across all groups simultaneously}
\usage{
test_all_groups_as(pbas, groupby, parallel = FALSE)
}
\arguments{
\item{pbas}{SummarizedExperiment: each row is a segment and each column is a pseudobulk.
Must contain assays named \code{i} (inclusion) and \code{e} (exclusion).
Rows correspond to segments; columns correspond to samples/pseudobulks.}

\item{groupby}{Either:
\itemize{
\item Vector of length \code{nrow(x)} or \code{ncol(x)} (for \code{SummarizedExperiment}): treated as the grouping factor directly.
\item Character scalar or vector of column names in \code{x} (or in \code{colData(x)} if \code{x} is a \code{SummarizedExperiment}):
those columns are concatenated row-wise (using \code{paste()} or \code{interaction()}) to form a grouping factor.
See \code{\link{get_groupby_factor}} for details.
}}

\item{parallel}{Logical; if \code{TRUE}, fit per‐segment GLMs in parallel (requires a registered \code{plyr} backend). Default \code{FALSE}.}
}
\value{
Data.frame with one row per segment (rownames = \code{rownames(pbas)} (segment IDs)), containing columns:
\itemize{
\item overdispersion: Estimated dispersion from each segment’s GLM
\item group: Raw p‐value from the likelihood‐ratio test for the \code{group} term
\item group_fdr: Benjamini‐Hochberg adjusted FDR (across all segments)
\item low_state: Group label with lowest mean PSI (from \code{\link{get_dpsi}})
\item high_state: Group label with highest mean PSI
\item dpsi: Difference in mean PSI between \code{high_state} and \code{low_state}
}
}
\description{
For a \code{SummarizedExperiment} of splicing segments with assays \code{i} (inclusion counts) and \code{e} (exclusion counts),
this function tests whether PSI differs across multiple groups (e.g., cell types) using a quasi‐binomial GLM.
It returns a data.frame of p‐values, adjusted FDR, and delta‐PSI for each segment.
}
\details{
\enumerate{
\item Converts \code{groupby} into a single grouping vector \code{group_factor} via \code{get_groupby_factor()}.
\item Calls \code{fit_as_glm()} with formula \code{x ~ group}, where \code{x} is the per‐segment \code{cbind(i,e)}.
\item Extracts raw p‐values for the \code{group} term, adjusts them (Benjamini-Hochberg) into \code{group_fdr}.
\item Computes delta‐PSI (\code{dpsi}), \code{low_state}, and \code{high_state} for each segment via \code{\link{get_dpsi}}.
\item Combine results (endure matching row order)
}
}
\seealso{
\code{\link{fit_as_glm}}, \code{\link{get_groupby_factor}}, \code{\link{get_dpsi}}, \code{\link{test_pair_as}}
}
