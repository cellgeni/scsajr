options(error = function(e) quit("no", 1))

delimeter <- "$"


#' Filter segments and samples based on coverage and group criteria
#'
#' This function takes a \code{SummarizedExperiment} containing alternative splicing data
#' (with assays “i” for inclusion counts and “e” for exclusion counts), and filters out:
#'   1. pseudobulk samples with fewer than \code{sample_min_ncells} cells,
#'   2. cell‐type groups with fewer than \code{celltype_min_samples} samples,
#'   3. segments (rows) that do not meet minimum coverage across pseudobulks,
#'   4. segments with low variability (standard deviation of PSI).
#'
#' @param pbas A \code{SummarizedExperiment} object. It must contain assays named “i”, “e”, and “psi”.
#'   Additionally, \code{colData(pbas)} should have a column (or combination of columns) specified by \code{groupby}, and a column \code{ncells} indicating how many cells contributed to each pseudobulk.
#' @param sample_filter Logical vector (length = number of columns in \code{pbas}) indicating which samples to retain initially.
#'   Defaults to \code{TRUE} for all samples.
#' @param sites Character vector of splice‐site patterns to keep (e.g., \code{c("ad", "aa", "dd")}). Only segments whose \code{rowData(pbas)$sites} is in this set will be retained.
#' @param min_cov Minimum total junction coverage (\code{i + e}) for a pseudobulk to count toward segment inclusion. Default 10.
#' @param sample_min_ncells Minimum number of cells per pseudobulk; pseudobulk samples with fewer cells are filtered out. Default 20.
#' @param celltype_min_samples Minimum number of pseudobulk samples per group (as defined by \code{groupby}) to keep that group. Default 2.
#' @param seg_min_samples The absolute minimum number of pseudobulks (after filtering by \code{min_cov}) that a segment must be observed in. Default 4.
#' @param seg_min_samples_fraq Minimum fraction of total pseudobulk samples that must meet \code{min_cov} for a segment to be retained. Default 0.2.
#'   The effective minimum samples per segment is \code{max(seg_min_samples, ncol(pbas) * seg_min_samples_fraq)}.
#' @param seg_min_sd Minimum standard deviation of PSI (percent spliced‐in) across pseudobulks for a segment to be retained. Default 0.1.
#' @param groupby Either (a) a character vector of length \code{ncol(pbas)} giving an existing grouping factor for each column, or
#'   (b) a character scalar or vector of column names in \code{colData(pbas)} whose pasted‐together values define the grouping factor.
#'   See \code{\link{get_groupby_factor}} for details.
#'
#' @return A filtered \code{SummarizedExperiment}, containing only those pseudobulks and segments that satisfy all criteria.
#'   The returned object will always have assays “i”, “e” converted to dense matrices, and will include a “psi” assay (percent spliced‐in) for each retained segment/samples.
#'   The \code{rowData} will also contain a column “nna” (number of pseudobulks with \code{i+e >= min_cov}) and “sd” (standard deviation of PSI).
#'
#' @examples
#' \dontrun{
#' data(pbas)
#' filtered_se <- filter_segments_and_samples(
#'   pbas = pbas,
#'   sample_filter = TRUE,
#'   sites = c("ad", "aa", "dd"),
#'   min_cov = 10,
#'   sample_min_ncells = 20,
#'   celltype_min_samples = 2,
#'   seg_min_samples = 4,
#'   seg_min_samples_fraq = 0.2,
#'   seg_min_sd = 0.1,
#'   groupby = "celltype"
#' )
#' head(filtered_se)
#' }
#' @export
filter_segments_and_samples <- function(
    pbas,
    sample_filter = TRUE,
    sites = c("ad", "aa", "dd"),
    min_cov = 10,
    sample_min_ncells = 20,
    celltype_min_samples = 2,
    seg_min_samples = 4,
    seg_min_samples_fraq = 0.2,
    seg_min_sd = 0.1,
    groupby = "celltype") {
  ## 1. Filter out pseudobulk samples with fewer than sample_min_ncells cells
  # Fleter pseudobulk samples with fewer cells
  sample_filter <- sample_filter & (SummarizedExperiment::colData(pbas)$ncells >= sample_min_ncells)

  # Determine how many samples remain in each group (celltype)
  group_factor <- get_groupby_factor(pbas, groupby) # returns a factor vector (length = ncol)
  ct_counts <- table(group_factor[sample_filter]) # Count samples per group (only among those passing sample_filter)

  # Only keep samples whose group has at least 'celltype_min_samples' pseudobulks
  sample_filter <- sample_filter & (group_factor %in% names(ct_counts)[ct_counts >= celltype_min_samples])

  pbas <- pbas[, sample_filter] # Subset the SummarizedExperiment columns (samples)


  ## 2. Filter out cell‐type groups with fewer than celltype_min_samples} samples
  # Keep only rows where rowData(pbas)$sites is in the specified 'sites' vector
  seg_sel <- SummarizedExperiment::rowData(pbas)$sites %in% sites
  pbas <- pbas[seg_sel, ]

  # Compute, for each segment, how many pseudobulks have (i + e) >= min_cov
  #   - assay(pbas, "i") and assay(pbas, "e") give inclusion/exclusion counts
  #   - Convert these to dense matrices to simplify row/column sums
  SummarizedExperiment::assay(pbas, "i") <- as.matrix(SummarizedExperiment::assay(pbas, "i"))
  SummarizedExperiment::assay(pbas, "e") <- as.matrix(SummarizedExperiment::assay(pbas, "e"))

  # 'nna' = number of pseudobulks with total coverage >= min_cov
  SummarizedExperiment::rowData(pbas)$nna <- rowSums((SummarizedExperiment::assay(pbas, "i") + SummarizedExperiment::assay(pbas, "e")) >= min_cov)


  ## 3. Filter out segments (rows) that do not meet minimum coverage across pseudobulks
  # Determine per‐segment minimum required pseudobulks:
  # If fractional rule (seg_min_samples_fraq) yields more required samples than seg_min_samples,
  # use the larger of the two.
  min_samples_required <- max(seg_min_samples, ceiling(ncol(pbas) * seg_min_samples_fraq))
  seg_sel2 <- SummarizedExperiment::rowData(pbas)$nna >= min_samples_required
  pbas <- pbas[seg_sel2, ]


  ## 4. Compute PSI (percent spliced‐in) and filter by segment variability
  # Calculate PSI = i / (i + e), setting PSI to NA whenever total coverage < min_cov
  psi_mat <- calc_psi(pbas, min_cov = min_cov)
  SummarizedExperiment::assay(pbas, "psi") <- psi_mat

  # Compute standard deviation of PSI across pseudobulks for each segment
  SummarizedExperiment::rowData(pbas)$sd <- apply(SummarizedExperiment::assay(pbas, "psi"), 1, stats::sd, na.rm = TRUE)

  # Keep only segments with standard deviation above threshold
  seg_sel3 <- SummarizedExperiment::rowData(pbas)$sd > seg_min_sd
  pbas <- pbas[seg_sel3, ]

  # Return the filtered SummarizedExperiment
  return(pbas)
}


#' Select marker segments per cell type
#'
#' Filters markers by false discovery rate (FDR) and absolute delta PSI (dPSI),
#' then selects the top n markers per cell type, and removing duplicates
#' (keeping only the strongest marker per segment across all cell_type).
#'
#' @param segment_stats A list containing numeric matrices `pv`, `fdr`, and `dpsi`
#'                      with rows = segments and cols = cell types.
#' @param n Integer: maximum number of markers to return per cell type (default: 5).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#' @param clean_duplicates Logical: if TRUE, keep only one entry per segment
#'                         (the one with highest |dpsi|) across all cell types (default: TRUE).
#' @return A data.frame with columns `pv`, `fdr`, `dpsi`, `seg_id`, and `group`(cell type),
#'         sorted by decreasing |dpsi|.
#' @examples
#' # Load example data
#' data(pbasf)
#'
#' # Select top 10 markers with stricter thresholds:
#' df <- select_markers(pbasf@metadata$markers, n = 10, fdr_thr = 0.01, dpsi_thr = 0.2)
#' head(df)
#' @export
select_markers <- function(segment_stats, n = 5, fdr_thr = 0.05, dpsi_thr = 0.1, clean_duplicates = TRUE) {
  # dataframe to store selected markers
  selected <- NULL

  # Loop over each cell_type (column in the FDR matrix)
  for (cell_type in colnames(segment_stats$fdr)) {
    # Filter segments based on FDR and dPSI thresholds
    is_significant <-
      (segment_stats$fdr[, cell_type] <= fdr_thr) &
        (abs(segment_stats$dpsi[, cell_type]) >= dpsi_thr)
    is_significant[is.na(is_significant)] <- FALSE # Treat NA values as non-significant

    # Skip if no segments pass for this cell_type
    if (!any(is_significant)) next

    # Build dataframe of all significant segments for this cell_type
    top <- data.frame(
      pv = segment_stats$pv[is_significant, cell_type],
      fdr = segment_stats$fdr[is_significant, cell_type],
      dpsi = segment_stats$dpsi[is_significant, cell_type],
      seg_id = rownames(segment_stats$pv)[is_significant],
      group = cell_type
    )

    # Order by decending |dPSI| and select top n segments
    keep_idx <- order(abs(top$dpsi), decreasing = TRUE)[seq_len(min(n, nrow(top)))]
    top <- top[keep_idx, ]

    # Append to output dataframe
    selected <- rbind(selected, top)
  }

  # Reset row name
  rownames(selected) <- NULL

  # (Optional) remove duplicate segment IDs (keep the highest |dpsi|)
  if (clean_duplicates) {
    selected <- lapply(split(selected, selected$seg_id), function(x) x[order(abs(x$dpsi), decreasing = TRUE)[1], ])
    selected <- do.call(rbind, selected)
  }

  # Final sort by descending |dpsi|
  selected <- selected[order(abs(selected$dpsi), decreasing = TRUE), ]

  # Return dataframe
  selected
}
