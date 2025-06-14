% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scsajr.R
\name{pseudobulk_metadata}
\alias{pseudobulk_metadata}
\title{Aggregate column‐level metadata for pseudobulk groups}
\usage{
pseudobulk_metadata(meta, groupby, aggregate = list(ncells = sum))
}
\arguments{
\item{meta}{A data.frame whose rows correspond to individual cells or samples,
and whose columns contain metadata (e.g., \code{sample_id}, \code{celltype}, \code{ncells}, etc.).
Row names should match column names of the corresponding \code{SummarizedExperiment} if applicable.}

\item{groupby}{Either:
\itemize{
\item Vector of length \code{nrow(x)} or \code{ncol(x)} (for \code{SummarizedExperiment}): treated as the grouping factor directly.
\item Character scalar or vector of column names in \code{x} (or in \code{colData(x)} if \code{x} is a \code{SummarizedExperiment}):
those columns are concatenated row-wise (using \code{paste()} or \code{interaction()}) to form a grouping factor.
See \code{\link{get_groupby_factor}} for details.
}}

\item{aggregate}{A named list of functions to apply to each column within each group.
For example, \code{list(ncells = sum, ngenes = max)} would sum the \code{ncells} column and take
the maximum of the \code{ngenes} column for each group. Default is \code{list(ncells = sum)}.}
}
\value{
A data.frame with one row per unique group. Columns include:
\itemize{
\item Any original columns in \code{meta} that had the same (constant) value for all rows in that group.
\item For each name in \code{aggregate}, a new column where the corresponding function has been applied
to that column’s values within each group.
Row names of the returned data.frame are the group labels.
}
}
\description{
Given a data.frame of per‐cell metadata and a grouping vector, this function aggregates metadata so that each row corresponds to one group.
For each group:
\enumerate{
\item Columns that are constant within the group are retained.
\item For any column listed in \code{aggregate}, a summary function is applied (e.g., sum).
}
}
\examples{
\dontrun{
# Suppose 'cell_meta' is a data.frame of 1000 cells with columns:
#   sample_id, celltype, ncells, batch, library_size
# We want to pseudobulk by 'celltype', summing 'ncells' and taking max 'library_size':
group_meta <- pseudobulk_metadata(
  meta = cell_meta,
  groupby = "celltype",
  aggregate = list(ncells = sum, library_size = max)
)
# 'group_meta' now has one row per unique celltype.
}

}
\seealso{
\code{\link{pseudobulk}}, \code{\link{get_groupby_factor}}
}
