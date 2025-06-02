#' PBASF: RangedSummarizedExperiment from 10X.GBM data
#'
#' A `RangedSummarizedExperiment` object generated from the 10X.GBM dataset
#' (processed via the cellgeni/nf-scsajr pipeline). It contains genomic ranges,
#' assay matrices, and associated metadata used throughout the package.
#'
#' @format A `RangedSummarizedExperiment` with 6 slots, including:
#' \describe{
#'   \item{rowRanges}{A `GRanges` of 6,794 features (exons and splice sites).}
#'   \item{colData}{A `DataFrame` of 10 samples with `sample_id`, `celltype`, etc.}
#'   \item{assays}{A list of matrices (`i`, `e`, `psi`) for counts and PSI values.}
#'   \item{metadata}{A list containing sub‐objects:
#'     \itemize{
#'       \item \code{markers}: marker tables per cell type
#'       \item \code{all_celltype_test}: differential testing results
#'       \item \code{go}, \code{ipro}: enrichment results
#'     }
#'   }
#' }
#' @source 10X.GBM dataset processed via \href{https://github.com/cellgeni/nf-scsajr}{cellgeni/nf-scsajr}
#' @keywords datasets
#' @name pbasf
#' @docType data
NULL
