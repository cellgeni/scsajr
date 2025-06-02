options(error = function(e) quit("no", 1))

delimeter <- "$"

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
