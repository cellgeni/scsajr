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
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
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


#' Fit a quasi‐binomial GLM for alternative splicing per segment
#'
#' For each segment (row) in a \code{SummarizedExperiment} with assays “i” (inclusion counts) and “e” (exclusion counts),
#' this function fits a quasi‐binomial generalized linear model of the form \code{i / (i + e) ~ predictors},
#' adding a pseudocount if desired. It can optionally return p‐values via a likelihood‐ratio test.
#'
#' @param pbas A \code{SummarizedExperiment} where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param formula An \code{\link[stats]{formula}} object specifying the model, e.g., \code{i ~ group}.
#'   Internally, the function treats “i” and “e” as a two‐column response.
#' @param data_terms A named list of covariates: data.frame columns referenced in \code{formula}.
#'   For example, if \code{formula = i ~ group}, then \code{data_terms = list(group = group_factor)} where \code{group_factor} is a vector of length \code{ncol(pbas)}.
#' @param pseudocount A numeric fraction of each count total to add to both “i” and “e”, expressed as a fraction of the total counts. Default 0 (no pseudocount).
#' @param parallel Logical; if \code{TRUE}, uses \code{\link[plyr]{alply}(…, .parallel = TRUE)} to fit models in parallel (requires registered backend). Default \code{FALSE}.
#' @param progress Character string passed to \code{.progress} in \code{\link[plyr]{alply}}. Default \code{"none"}.
#' @param return_pv Logical; if \code{TRUE}, perform a quasi‐likelihood ratio test per segment to return p‐values for each term in \code{formula}.
#'   Otherwise, return fitted \code{\link[stats]{glm}} objects. Default \code{FALSE}.
#' @param overdisp Logical; if \code{TRUE}, account for overdispersion when computing test statistics. Default \code{TRUE}.
#' @param disp_param Optional numeric vector of length \code{nrow(pbas)} providing dispersion parameters per segment.
#'   If \code{NULL}, dispersion is estimated from each model’s residuals.
#'
#' @return If \code{return_pv = FALSE}, a named list of length \code{nrow(pbas)} of fitted \code{glm} objects
#'   (or \code{NA} for segments where fitting failed). Names correspond to \code{rownames(pbas)}.
#'   If \code{return_pv = TRUE}, a data.frame with one row per segment, columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion (or provided \code{disp_param}).}
#'     \item{<term>}{Likelihood‐ratio p‐value for each term in \code{formula}.}
#'   }
#'   Row names correspond to \code{rownames(pbas)}.
#'
#' @details
#' The function loops over each segment, constructs a two‐column matrix \code{cbind(i, e)}, adds pseudocounts if requested, and fits a quasi‐binomial \code{\link[stats]{glm}}.
#' If \code{return_pv = TRUE}, it calls \code{\link{asq_lrt}} to compute a likelihood‐ratio test for each term in \code{formula}.
#' Any errors or warnings during fitting are caught; segments with fitting errors return \code{NA}.
#'
#' @seealso \code{\link{asq_lrt}}, \code{\link[stats]{glm}}, \code{\link[plyr]{alply}}
#' @export
fit_as_glm <- function(
    pbas,
    formula,
    data_terms,
    pseudocount = 0,
    parallel = FALSE,
    progress = "none",
    return_pv = FALSE,
    overdisp = TRUE,
    disp_param = NULL) {
  # Ensure plyr is available
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Please install the 'plyr' package to use fit_as_glm().")
  }

  # Extract term names from the formula
  term_labels <- attr(terms(formula), "term.labels")

  # For each segment (row), fit a quasi-binomial GLM
  results_list <- plyr::alply(
    .data = seq_len(nrow(pbas)),
    .margins = 1,
    .fun = function(i) {
      # 1. Build the i/e matrix and add pseudocounts if specified
      inc_vec <- SummarizedExperiment::assay(pbas, "i")[i, ]
      exc_vec <- SummarizedExperiment::assay(pbas, "e")[i, ]
      tot <- inc_vec + exc_vec

      # Add pseudocount proportional to total counts for each sample
      if (pseudocount > 0) {
        adj <- tot * pseudocount
        inc_vec <- inc_vec + adj
        exc_vec <- exc_vec + adj
      }

      # Place adjusted counts into the data frame used by glm()
      #    'x' must be a two-column matrix: cbind(i, e)
      data_terms$x <- cbind(inc_vec, exc_vec)

      # 2. Fit the quasi-binomial GLM
      fit <- tryCatch(
        stats::glm(formula, data = data_terms, family = "quasibinomial"),
        error = function(e) {
          warning(paste0("Segment ", rownames(pbas)[i], ": ", e$message))
          return(NA)
        },
        warning = function(w) {
          warning(paste0("Segment ", rownames(pbas)[i], ": ", w$message))
          return(NA)
        }
      )

      # 3. If return_pv and fit succeeded, compute p-values via LRT
      if (return_pv && !is.na(fit)[1]) {
        fit <- asq_lrt(
          fit,
          overdisp = overdisp,
          disp_param = disp_param[i],
          term_labels = term_labels,
          seg_id = rownames(pbas)[i]
        )
      }

      fit
    },
    .parallel = parallel,
    .progress = progress
  )

  # If return_pv, combine p-value results into a data.frame
  if (return_pv) {
    pv_df <- do.call(rbind, results_list)
    rownames(pv_df) <- rownames(pbas)
    colnames(pv_df) <- c("overdispersion", term_labels)
    return(pv_df)
  }

  # Otherwise, return a named list of glm objects (or NA)
  names(results_list) <- rownames(pbas)
  attr(results_list, "term.labels") <- term_labels
  return(results_list)
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
