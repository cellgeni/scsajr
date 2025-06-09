#!/usr/bin/env Rscript

# Need to add package description
# Need to figureout how to keep links in the Roxygen documentation

#' Filter segments and samples based on coverage and group criteria
#'
#' This function takes a `SummarizedExperiment` containing alternative splicing data
#' (with assays “i” for inclusion counts and “e” for exclusion counts), and filters out:
#'   1. pseudobulk samples with fewer than `sample_min_ncells` cells,
#'   2. cell‐type groups with fewer than `celltype_min_samples` samples,
#'   3. segments (rows) that do not meet minimum coverage across pseudobulks,
#'   4. segments with low variability (standard deviation of PSI).
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param sample_filter Logical vector (length = number of columns in `pbas`) indicating which samples to retain initially.
#'   Defaults to `TRUE` for all samples.
#' @param sites Character vector of splice‐site patterns to keep (e.g., `c("ad", "aa", "dd")`). Only segments whose `rowData(pbas)$sites` is in this set will be retained.
#' @param min_cov Minimum total junction coverage (`i + e`) for a pseudobulk to count toward segment inclusion. Default 10.
#' @param sample_min_ncells Minimum number of cells per pseudobulk; pseudobulk samples with fewer cells are filtered out. Default 20.
#' @param celltype_min_samples Minimum number of pseudobulk samples per group (as defined by `groupby`) to keep that group. Default 2.
#' @param seg_min_samples The absolute minimum number of pseudobulks (after filtering by `min_cov`) that a segment must be observed in. Default 4.
#' @param seg_min_samples_fraq Minimum fraction of total pseudobulk samples that must meet `min_cov` for a segment to be retained. Default 0.2.
#'   The effective minimum samples per segment is `max(seg_min_samples, ncol(pbas) * seg_min_samples_fraq)`.
#' @param seg_min_sd Minimum standard deviation of PSI (percent spliced‐in) across pseudobulks for a segment to be retained. Default 0.1.
#' @param groupby Either (a) a character vector of length `ncol(pbas)` giving an existing grouping factor for each column, or
#'   (b) a character scalar or vector of column names in `colData(pbas)` whose pasted‐together values define the grouping factor.
#'   See `\link{get_groupby_factor}` for details.
#'
#' @return A filtered `SummarizedExperiment`, containing only those pseudobulks and segments that satisfy all criteria.
#'   The returned object will always have assays “i”, “e” converted to dense matrices, and will include a “psi” assay (percent spliced‐in) for each retained segment/samples.
#'   The `rowData` will also contain a column “nna” (number of pseudobulks with `i+e >= min_cov`) and “sd” (standard deviation of PSI).
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
  # Filter pseudobulk samples with fewer cells
  sample_filter <- sample_filter & (SummarizedExperiment::colData(pbas)$ncells >= sample_min_ncells)

  # Determine how many samples remain in each group (celltype)
  group_factor <- get_groupby_factor(pbas, groupby) # returns a factor vector (length = ncol)
  ct_counts <- table(group_factor[sample_filter]) # Count samples per group (only among those passing sample_filter)

  # Only keep samples whose group has at least 'celltype_min_samples' pseudobulks
  sample_filter <- sample_filter & (group_factor %in% names(ct_counts)[ct_counts >= celltype_min_samples])

  pbas <- pbas[, sample_filter] # Subset the SummarizedExperiment columns (samples)


  ## 2. Filter out cell‐type groups with fewer than celltype_min_samples samples
  # Keep only rows where rowData(pbas)$sites is in the specified 'sites' vector
  seg_sel <- SummarizedExperiment::rowData(pbas)$sites %in% sites
  pbas <- pbas[seg_sel, ]

  # Compute, for each segment, how many pseudobulks have (i + e) >= min_cov
  #   - assay(pbas, "i") and assay(pbas, "e") give inclusion/exclusion counts
  #   - Convert these to dense matrices to simplify row/column sums
  SummarizedExperiment::assay(pbas, "i") <- as.matrix(SummarizedExperiment::assay(pbas, "i"))
  SummarizedExperiment::assay(pbas, "e") <- as.matrix(SummarizedExperiment::assay(pbas, "e"))

  # 'nna' = number of pseudobulks with total coverage >= min_cov
  SummarizedExperiment::rowData(pbas)$nna <- Matrix::rowSums((SummarizedExperiment::assay(pbas, "i") + SummarizedExperiment::assay(pbas, "e")) >= min_cov)


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
#' For each segment (row) in a `SummarizedExperiment` with assays “i” (inclusion counts) and “e” (exclusion counts),
#' this function fits a quasi‐binomial generalized linear model of the form `i / (i + e) ~ predictors`,
#' adding a pseudocount if desired. It can optionally return p‐values via a likelihood‐ratio test.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param formula An `\link[stats]{formula}` object specifying the model, e.g., `i ~ group`.
#'   Internally, the function treats “i” and “e” as a two‐column response.
#' @param data_terms A named list of covariates: data.frame columns referenced in `formula`.
#'   For example, if `formula = i ~ group`, then `data_terms = list(group = group_factor)` where `group_factor` is a vector of length `ncol(pbas)`.
#' @param pseudocount A numeric fraction of each count total to add to both “i” and “e”, expressed as a fraction of the total counts. Default 0 (no pseudocount).
#' @param parallel Logical; if `TRUE`, uses `\link[plyr]{alply}(…, .parallel = TRUE)` to fit models in parallel (requires registered backend). Default `FALSE`.
#' @param progress Character string passed to `.progress` in `\link[plyr]{alply}`. Default `"none"`.
#' @param return_pv Logical; if `TRUE`, perform a quasi‐likelihood ratio test per segment to return p‐values for each term in `formula`.
#'   Otherwise, return fitted `\link[stats]{glm}` objects. Default `FALSE`.
#' @param overdisp Logical; if `TRUE`, account for overdispersion when computing test statistics. Default `TRUE`.
#' @param disp_param Optional numeric vector of length `nrow(pbas)` providing dispersion parameters per segment.
#'   If `NULL`, dispersion is estimated from each model’s residuals.
#'
#' @return If `return_pv = FALSE`, a named list of length `nrow(pbas)` of fitted `glm` objects
#'   (or `NA` for segments where fitting failed). Names correspond to `rownames(pbas)`.
#'   If `return_pv = TRUE`, a data.frame with one row per segment, columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion (or provided `disp_param`).}
#'     \item{<term>}{Likelihood‐ratio p‐value for each term in `formula`.}
#'   }
#'   Row names correspond to `rownames(pbas)`.
#'
#' @details
#' The function loops over each segment, constructs a two‐column matrix `cbind(i, e)`, adds pseudocounts if requested, and fits a quasi‐binomial `\link[stats]{glm}`.
#' If `return_pv = TRUE`, it calls `\link{qbinom_lrt}` to compute a likelihood‐ratio test for each term in `formula`.
#' Any errors or warnings during fitting are caught; segments with fitting errors return `NA`.
#'
#' @seealso `\link{qbinom_lrt}`, `\link[stats]{glm}`, `\link[plyr]{alply}`
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
  # Extract term names from the formula
  term_labels <- attr(stats::terms(formula), "term.labels")

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
        fit <- qbinom_lrt(
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


#' Quasi‐likelihood ratio test (Chi-square ANOVA) for a fitted quasi‐binomial GLM
#'
#' Given a fitted `glm(..., family = "quasibinomial")`, this function computes a likelihood‐ratio test for each term in the model, accounting for overdispersion.
#'
#' @param fit_obj A fitted `glm` object (quasi‐binomial family).
#' @param overdisp Logical; if `TRUE`, use estimated dispersion from `fit_obj` (or `disp_param`) when computing test statistics.
#'   If `FALSE`, force dispersion = 1 (i.e., treat as simple binomial).
#' @param disp_param Optional numeric; if provided, use this as the dispersion parameter for the test.
#'   Otherwise, estimate dispersion from `fit_obj`’s residuals.
#' @param term_labels Character vector of term names (excluding the intercept) in the original formula, in the same order returned by `attr(terms(formula), "term.labels")`.
#' @param seg_id Character; segment identifier (row name) used in warning messages.
#'
#' @return A numeric vector of length `length(term_labels) + 1`.
#'   The first element is the dispersion estimate actually used (either `disp_param` or `summary(fit_obj)$dispersion`),
#'    and the remaining elements are p‐values (Chisq tests) for each term in `term_labels`, in the same order.
#'   If `fit_obj` is `NA` (fitting failed), returns a vector of `NA`s of appropriate length.
#'
#' @details
#' 1. If `disp_param` is `NULL`, we take `summary(fit_obj)$dispersion` as the overdispersion estimate.
#'    If the model has zero residual degrees of freedom, dispersion cannot be estimated; in that case, we return `NA` for dispersion and issue a warning (unless `overdisp = FALSE`).
#' 2. We set `disp_used = max(1, <disp_estimate>)` so that the effective dispersion is at least 1. If `overdisp = FALSE`, we force `disp_used = 1`.
#' 3. We then call `\link[stats]{anova}(fit_obj, test = "Chisq", dispersion = disp_used)`.
#'    The p‐values for each term appear in column 5 of the ANOVA table.
#' 4. If any error occurs when computing the ANOVA, we return `rep(NA, length(term_labels) + 1)`.
#'
#' @seealso `\link{fit_as_glm}`, `\link[stats]{anova}`, `\link[stats]{glm}`
#' @export
qbinom_lrt <- function(
    fit_obj,
    overdisp = TRUE,
    disp_param = NULL,
    term_labels,
    seg_id) {
  # Number of output elements: 1 for dispersion, plus one per term
  n_out <- length(term_labels) + 1
  result <- rep(NA_real_, n_out)

  # If fitting failed (NA), return all NAs
  if (is.na(fit_obj)[1]) {
    return(result)
  }

  # Determine dispersion estimate
  if (is.null(disp_param)) {
    disp_est <- summary(fit_obj)$dispersion
    # If no residual DF, cannot estimate dispersion
    if (fit_obj$df.residual == 0) {
      disp_est <- NA_real_
      if (overdisp) {
        warning(paste0(
          "Segment ", seg_id,
          ": cannot estimate overdispersion without replicates; returning NA"
        ))
      }
    }
  } else {
    disp_est <- disp_param
  }

  # Use at least 1 for dispersion; if overdisp = FALSE, force 1
  disp_used <- max(1, disp_est, na.rm = TRUE)
  if (!overdisp) {
    disp_used <- 1
  }
  result[1] <- disp_used

  # Run ANOVA (Chisq) with specified dispersion; catch errors
  a_tbl <- tryCatch(
    stats::anova(fit_obj, test = "Chisq", dispersion = disp_used),
    error = function(e) {
      warning(paste0(
        "Segment ", seg_id,
        ": ANOVA Chi-square test failed: ", e$message
      ))
      return(NULL)
    }
  )
  if (is.null(a_tbl)) {
    return(result)
  }

  # Extract p‐values for each term: they occupy rows 2:(1+length(term_labels)), column 5
  # (anova table has rows: intercept, term1, term2, ..., Residuals)
  pvals <- a_tbl[2:(1 + length(term_labels)), 5]
  result[-1] <- pvals
  return(result)
}


#' Determine grouping factor from a data.frame or SummarizedExperiment
#'
#' If `groupby` is a vector (e.g., a factor) whose length matches the number of rows of a `data.frame` or the number of columns of a `SummarizedExperiment`, it is returned unchanged.
#' Otherwise, if `groupby` is a character scalar or vector of column names present in the input,
#'  those columns are pasted together row-wise (for a `data.frame`) or via `interaction()` (for factor columns in`SummarizedExperiment::colData()`) to form a single grouping factor.
#'
#' @param x Either a `SummarizedExperiment` (in which case its `colData()` is used) or a `data.frame` whose rows correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `nrow(x)` or `ncol(x)` (for `SummarizedExperiment`):: treated as the grouping factor directly.
#'     \item A character scalar or vector of column names in `x` (or in `colData(x)` if `x` is a `SummarizedExperiment`):
#'            those columns are concatenated row-wise (using `paste()` or `interaction()`) to form a grouping factor.
#'   }
#' @param sep Character string used to separate pasted values when `groupby` is multiple column names. Default is `"\$"`.
#'
#' @return A character vector (or factor) of length `nrow(x)`, representing the grouping factor.
#'   If `groupby` was already the same length as `nrow(x)`, it is returned unchanged.
#'   Otherwise, it pastes together columns of `x`.
#'
#' @return A vector (typically a factor or character vector) of length equal to `nrow(x)` (for `data.frame`) or `ncol(x)` (for `SummarizedExperiment`), representing the grouping factor.
#'         If `groupby` was already the correct length, it is returned unchanged.
#'         Otherwise, it pastes or interacts specified columns of `x`.
#'
#' @examples
#' df <- data.frame(
#'   sample = paste0("S", 1:4),
#'   celltype = c("A", "A", "B", "B"),
#'   batch = c("X", "Y", "X", "Y")
#' )
#'
#' # Case 1: groupby is a factor vector of correct length
#' grp1 <- get_groupby_factor(df, groupby = c("A", "A", "B", "B"))
#'
#' # Case 2: groupby is one column name
#' grp2 <- get_groupby_factor(df, groupby = "celltype")
#'
#' # Case 3: groupby is multiple column names
#' grp3 <- get_groupby_factor(df, groupby = c("celltype", "batch"))
#'
#' # When x is a SummarizedExperiment
#' # sce$group_col <- factor(paste0(sce$sample_id, "_", sce$celltype))
#' # grp_se <- get_groupby_factor(sce, groupby = "group_col")
#'
#' print(grp1)
#' print(grp2)
#' print(grp3)
#'
#' @seealso `\link{filter_segments_and_samples}`, `\link{test_all_groups_as}`
#' @export
get_groupby_factor <- function(x, groupby, sep = "$") {
  # Determine whether x is a SummarizedExperiment or a data.frame
  is_se <- inherits(x, "SummarizedExperiment")

  # If SummarizedExperiment, extract its colData
  if (is_se) {
    coldat <- SummarizedExperiment::colData(x)
    n <- ncol(x) # Number of samples / columns
  } else if (is.data.frame(x)) {
    coldat <- x
    n <- nrow(x) # Number of samples / rows
  } else {
    stop(
      "get_groupby_factor(): `x` must be either a data.frame or a SummarizedExperiment."
    )
  }

  # 1. If groupby is a vector/factor of correct length, return it unchanged
  if ((is.atomic(groupby) || is.factor(groupby)) && length(groupby) == n) {
    return(groupby)
  }

  # 2. If groupby is character vector of column names in coldat/data.frame
  if (is.character(groupby)) {
    # (a) If all elements of groupby match column names in coldat
    if (all(groupby %in% colnames(coldat))) {
      # Extract requested columns
      subset_df <- as.data.frame(coldat[, groupby, drop = FALSE])

      # If only one column, return that column directly (as factor or character)
      if (length(groupby) == 1) {
        return(subset_df[[1]])
      }

      # More than one column: if all are factors, use interaction(); else, paste()
      are_factors <- vapply(subset_df, is.factor, logical(1))
      if (all(are_factors)) {
        # Combine factor columns into one factor via interaction
        combined_factor <- interaction(subset_df, drop = TRUE, sep = sep)
        return(combined_factor)
      } else {
        # Convert all columns to character and paste row-wise
        pasted <- do.call(paste, c(lapply(subset_df, as.character), sep = sep))
        return(pasted)
      }
    }
  }

  # If we reach here, groupby was neither a length-n vector nor valid column names
  stop(
    "`groupby` must be either:\n",
    " - A vector/factor of length ", n, ",\n",
    " - A character scalar or vector of column names in the data."
  )
}


#' Test alternative splicing across all groups simultaneously
#'
#' For a `SummarizedExperiment` of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests whether PSI differs across multiple groups (e.g., cell types) using a quasi‐binomial GLM.
#' It returns a data.frame of p‐values, adjusted FDR, and delta‐PSI for each segment.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segment IDs; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#'
#' @return A data.frame with one row per segment (rownames = `rownames(pbas)` (segment IDs)), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a single grouping vector `group_factor` via `get_groupby_factor()`.
#' 2. Calls `fit_as_glm()` with formula `x ~ group`, where `x` is the per‐segment `cbind(i,e)`.
#' 3. Extracts raw p‐values for the `group` term, adjusts them (Benjamini-Hochberg) into `group_fdr`.
#' 4. Computes delta‐PSI (`dpsi`), `low_state`, and `high_state` for each segment via `get_dpsi()`.
#' 5. Combine results (endure matching row order)
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, `\link{get_dpsi}`, `\link{test_pair_as}`
#' @export
test_all_groups_as <- function(
    pbas,
    groupby,
    parallel = FALSE) {
  # 1. Construct grouping factor
  group_factor <- get_groupby_factor(pbas, groupby)

  # 2. Fit GLMs per segment, returning a data.frame of (overdispersion + p‐values)
  pv_df <- fit_as_glm(
    pbas = pbas,
    formula = x ~ group,
    data_terms = list(group = group_factor),
    return_pv = TRUE,
    parallel = parallel
  )
  # pv_df has columns: overdispersion, group (p‐value)

  # 3. Adjust p‐values (Benjamini-Hochberg)
  pv_df <- as.data.frame(pv_df, stringsAsFactors = FALSE)
  pv_df$group_fdr <- stats::p.adjust(pv_df$group, method = "BH")

  # 4. Compute delta‐PSI, low_state, high_state
  dpsi_df <- get_dpsi(pbas, group_factor, min_cov = 50)

  # 5. Combine results: ensure matching row order
  res <- cbind(
    overdispersion = pv_df$overdispersion,
    group = pv_df$group,
    group_fdr = pv_df$group_fdr,
    dpsi_df[c("low_state", "high_state", "dpsi")]
  )
  rownames(res) <- rownames(pbas)
  return(res)
}


#' Test alternative splicing between two conditions (pairwise)
#'
#' For a `SummarizedExperiment` of splicing segments with assays “i” (inclusion counts) and “e” (exclusion counts),
#'  this function tests for differential splicing between two specified groups (e.g., two cell types) using a quasi‐binomial GLM per segment.
#' It returns a data.frame of p‐values, FDR, and delta‐PSI between the two conditions.
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays named “i” (inclusion) and “e” (exclusion).
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param conditions A length‐2 character vector specifying the two group labels to compare (e.g., `c("A", "B")`).
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#'
#' @return A data.frame with one row per segment (rownames = `rownames(pbas)`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{pv}{Raw p‐value for the group term (two‐level factor).}
#'     \item{fdr}{Benjamini-Hochberg FDR (across all segments).}
#'     \item{dpsi}{Difference in mean PSI between `conditions[2]` and `conditions[1]`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a single grouping vector `group_factor` via `get_groupby_factor()`.
#' 2. Filters `group_factor` and `pbas` to keep only columns corresponding to the two `conditions`.
#' 3. Calls `fit_as_glm()` with `formula = x ~ group` to get raw p‐values for each segment.
#' 4. Adjusts p‐values (Benjamini-Hochberg) into `fdr`.
#' 5. Computes delta‐PSI (`dpsi`) per segment by pseudobulking and calling `calc_psi()`, then taking
#'    mean PSI per condition and subtracting.
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, \code{\link{calc_psi}}, \link{pseudobulk}, `\link{test_all_groups_as}`
#' @export
test_pair_as <- function(
    pbas,
    groupby,
    conditions,
    parallel = FALSE) {
  # 1. Construct grouping factor
  group_factor <- get_groupby_factor(pbas, groupby)

  # 2. Identify columns (samples) that belong to the two conditions
  keep_cols <- which(group_factor %in% conditions)
  group_factor <- group_factor[keep_cols]
  pbas_sub <- pbas[, keep_cols]

  # 3. Fit GLMs per segment on only the two conditions, returning raw p‐values
  pv_df <- fit_as_glm(
    pbas = pbas_sub,
    formula = x ~ group,
    data_terms = list(group = group_factor),
    return_pv = TRUE,
    parallel = parallel
  )
  # pv_df has columns: overdispersion, group

  # 4. Adjust p‐values (Benjamini-Hochberg)
  pv_df <- as.data.frame(pv_df, stringsAsFactors = FALSE)
  pv_df$fdr <- stats::p.adjust(pv_df$group, method = "BH")

  # 5. Compute PSI matrix for the two‐condition pseudobulk
  pb <- pseudobulk(pbas_sub, group_factor)
  psi_mat <- calc_psi(pb)
  # Extract PSI for the two conditions in the order given by 'conditions'
  psi_vals <- psi_mat[, conditions, drop = FALSE]
  dpsi_vec <- psi_vals[, 2] - psi_vals[, 1]

  # 6. Combine results into final data.frame
  res <- data.frame(
    overdispersion = pv_df$overdispersion,
    pv = pv_df$group,
    fdr = pv_df$fdr,
    dpsi = dpsi_vec,
    row.names = rownames(pbas_sub),
    stringsAsFactors = FALSE
  )

  # Ensure rows are in the same order as the original pbas
  res <- res[rownames(pbas), , drop = FALSE]
  return(res)
}


#' Identify marker segments for each group via one‐vs‐rest tests
#'
#' For a `SummarizedExperiment` of splicing segments with assays “psi” (percent spliced‐in),
#'  this function conducts, for each unique group label, a one‐vs‐rest quasi‐binomial GLM test
#'  (i.e., comparing that group to all other samples).
#' It returns p‐values and delta‐PSI matrices (segments × groups).
#'
#' @param pbas A `SummarizedExperiment` where each row is a segment and each column is a pseudobulk.
#'   Must contain assays “i”, “e”, and “psi”.
#'   Rows correspond to segments; columns correspond to samples/pseudobulks.
#' @param groupby Either:
#'   \itemize{
#'     \item A vector of length `ncol(pbas)` that directly gives a group label for each column, or
#'     \item One or more column names in `colData(pbas)`, whose values (pasted together if multiple) define the grouping factor.
#'   }
#'   See `\link{get_groupby_factor}` for details.
#' @param parallel Logical; if `TRUE`, fit per‐segment GLMs in parallel (requires a registered `plyr` backend). Default `FALSE`.
#' @param verbose Logical; if `TRUE`, prints each group being tested. Default `FALSE`.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{pv}{Matrix: Raw p‐values, with rows = segments, columns = group labels}
#'     \item{fdr}{Matrix: Benjamini-Hochberg FDR (across all segments)}
#'     \item{dpsi}{Matrix: Delta‐PSI (mean PSI(group) - mean PSI(others)) per segment}
#'   }
#'   Row names of each matrix are `rownames(pbas)`; column names are the unique groups.
#'
#' @details
#' 1. Builds a grouping factor via `get_groupby_factor(pbas, groupby)`.
#' 2. For each unique label `group` in that factor:
#'    \enumerate{
#'      \item Creates a binary factor `f` where `f == TRUE` if sample’s group == `group`, FALSE otherwise.
#'      \item Calls `fit_as_glm()` with `formula = x ~ f` and `data_terms = list(f = f)`, requesting p‐values (`return_pv = TRUE`).
#'      \item Extracts the “group” p‐value column (corresponding to `f`) for each segment and stores in `pv[, group]`.
#'      \item Computes `dpsi[, group]` as `mean(psi[segment, f]) - mean(psi[segment, !f])`.
#'    }
#' 3. Adjusts each column of `pv` by Benjamini-Hochberg into `fdr[, group]`.
#'
#' @seealso `\link{fit_as_glm}`, `\link{get_groupby_factor}`, `\link{calc_psi}`
#' @export
find_marker_as <- function(
    pbas,
    groupby,
    parallel = FALSE,
    verbose = FALSE) {
  # 1. Construct grouping vector (length = ncol)
  group_factor <- get_groupby_factor(pbas, groupby)
  unique_groups <- unique(group_factor)
  n_segs <- nrow(pbas)
  n_groups <- length(unique_groups)

  # 2. Initialise result matrices
  pv_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )
  dpsi_mat <- matrix(NA_real_,
    nrow = n_segs, ncol = n_groups,
    dimnames = list(rownames(pbas), unique_groups)
  )

  # 3. For each group, perform one‐vs‐rest test
  for (group in unique_groups) {
    if (verbose) {
      message("Testing group: ", group)
    }
    # Binary factor: TRUE if sample belongs to this group
    f <- group_factor == group

    # Fit GLM per segment; return a data.frame with columns: overdispersion and 'group' p‐value
    glm_res <- fit_as_glm(
      pbas       = pbas,
      formula    = x ~ f,
      data_terms = list(f = f),
      return_pv  = TRUE,
      parallel   = parallel
    )
    # glm_res[data.frame]: rownames = segment IDs, col 'group' = p‐value

    # Store raw p‐values for this group
    pv_mat[, group] <- glm_res[, "group"]

    # Compute dpsi: for each segment, mean PSI in 'group' minus mean PSI in others
    psi_vals <- SummarizedExperiment::assay(pbas, "psi")

    # For rows with all-NA in either subset, result is NA automatically
    mean_ct <- rowMeans(psi_vals[, f, drop = FALSE], na.rm = TRUE)
    mean_other <- rowMeans(psi_vals[, !f, drop = FALSE], na.rm = TRUE)
    dpsi_mat[, group] <- mean_ct - mean_other
  }

  # 4. Adjust p‐values column‐wise (BH)
  fdr_mat <- apply(pv_mat, 2, stats::p.adjust, method = "BH")

  # 5. Return list of matrices
  return(list(
    pv = pv_mat,
    fdr = fdr_mat,
    dpsi = dpsi_mat
  ))
}


#' Select top marker segments per group (celltype)
#'
#' From a list of per‐segment test results (as returned by `find_marker_as()`),
#'  this function selects up to `n` top segments per group
#'  whose FDR is below a threshold and whose absolute delta‐PSI exceeds a threshold.
#' Optionally, it removes duplicate segments so that each segment appears only once
#'  (assigned to the group where it has the largest |dpsi|).
#'
#' @param markers A list containing numeric matrices `pv`, `fdr`, and `dpsi`,
#'   each a matrix with rows = segments and columns = group labels (as returned by `find_marker_as()`).
#'   \describe{
#'     \item{`markers$pv`}{Matrix of raw p‐values, dimensions: ∣segments∣ × ∣groups∣.}
#'     \item{`markers$fdr`}{Matrix of Benjamini‐Hochberg adjusted FDR p‐values, same dimensions.}
#'     \item{`markers$dpsi`}{Matrix of delta‐PSI (mean PSI(group) - mean PSI(others)), same dimensions.}
#'   }
#' @param n Integer: maximum number of markers to return per cell type (default: 5).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#' @param clean_duplicates Logical: if `TRUE`, keep only one entry per segment
#'                         (the one with highest |dpsi|) across all cell types (default: `TRUE`).
#'
#' @return A data.frame with selected marker segments. Columns:
#'   \describe{
#'     \item{`pv`}{Raw p‐value from the likelihood‐ratio test for the group effect in the selected group.}
#'     \item{`fdr`}{Benjamini‐Hochberg adjusted FDR value for the segment in the selected group.}
#'     \item{`dpsi`}{Delta‐PSI (signed) for the segment in the selected group.}
#'     \item{`seg_id`}{Segment identifier (rownames of `markers$pv`).}
#'     \item{`group`}{Group label for which this segment is selected.}
#'   }
#'   Rows are ordered by descending |dpsi| across all selected segments.
#'
#' @details
#' For each group (column) in `markers$fdr`, the function:
#' \enumerate{
#'   \item Identifies segments `seg` satisfying
#'         `markers$fdr[seg, group] <= fdr_thr` and
#'         `|markers$dpsi[seg, group]| >= dpsi_thr`.
#'   \item If at least one segment passes, constructs a temporary data.frame `t` containing
#'         `pv = markers$pv[seg, group]`,
#'         `fdr = markers$fdr[seg, group]`,
#'         `dpsi = markers$dpsi[seg, group]`,
#'         `seg_id = seg`,
#'         `group = group`.
#'   \item Selects up to `n` rows from `t`, ordering by `abs(dpsi)` in descending order.
#'   \item Appends these rows to an aggregate `res` data.frame.
#' }
#' If `clean_duplicates = TRUE`, it then splits `res` by `seg_id` and retains only the row
#'  where `abs(dpsi)` is maximal across groups, removing duplicates.
#' Finally, it re‐orders `res` by `abs(dpsi)` in descending order and returns it.
#'
#' @examples
#' # Load example data
#' data(pbasf)
#'
#' # Select top 10 markers with stricter thresholds:
#' df <- select_markers(pbasf@metadata$markers, n = 10, fdr_thr = 0.01, dpsi_thr = 0.2)
#' head(df)
#'
#' @seealso `\link{find_marker_as}`, `\link{test_all_groups_as}`, `\link{test_pair_as}`
#' @export
select_markers <- function(
    markers,
    n = 5,
    fdr_thr = 0.05,
    dpsi_thr = 0.1,
    clean_duplicates = TRUE) {
  # Initialise an empty result
  result <- NULL

  # Iterate over each group (celltype) (column name) in the fdr matrix
  for (group in colnames(markers$fdr)) {
    # 1. Identify segments passing both FDR and |dpsi| thresholds
    pass_fdr <- markers$fdr[, group] <= fdr_thr
    pass_dpsi <- abs(markers$dpsi[, group]) >= dpsi_thr
    keep_idx <- which(pass_fdr & pass_dpsi)

    # Replace NA flags with FALSE
    keep_idx <- keep_idx[!is.na(keep_idx)]

    # If no segments pass, skip to next group
    if (length(keep_idx) == 0) {
      next
    }

    # 2. Build a temporary data.frame of candidate segments for this group
    temp_df <- data.frame(
      pv = markers$pv[keep_idx, group],
      fdr = markers$fdr[keep_idx, group],
      dpsi = markers$dpsi[keep_idx, group],
      seg_id = rownames(markers$pv)[keep_idx],
      group = group,
      stringsAsFactors = FALSE
    )

    # 3. Order by descending |dpsi| and take up to 'n' rows
    order_idx <- order(abs(temp_df$dpsi), decreasing = TRUE)
    top_n_idx <- utils::head(order_idx, n)
    temp_sel <- temp_df[top_n_idx, , drop = FALSE]

    # 4. Append to the cumulative result
    result <- if (is.null(result)) {
      temp_sel
    } else {
      rbind(result, temp_sel)
    }
  }

  # Reset row names of the aggregated result
  rownames(result) <- NULL

  # 5. If clean_duplicates, ensure each seg_id appears only once
  if (clean_duplicates && nrow(result) > 0) {
    # Split by seg_id; within each group, keep row with largest |dpsi|
    by_seg <- split(result, result$seg_id)
    cleaned_list <- lapply(by_seg, function(df_seg) {
      # Order rows for this segment by descending |dpsi|
      idx <- order(abs(df_seg$dpsi), decreasing = TRUE)
      df_seg[idx[1], , drop = FALSE]
    })
    result <- do.call(rbind, cleaned_list)
    rownames(result) <- NULL
  }

  # 6. Finally, order all rows by descending |dpsi|
  if (nrow(result) > 0) {
    result <- result[order(abs(result$dpsi), decreasing = TRUE), , drop = FALSE]
  }

  return(result)
}

#' Select markers from all‐celltype test results
#'
#' From a data.frame of per‐segment results comparing all groups simultaneously (e.g., output of `test_all_groups_as()`),
#'   this function selects only those segments that pass both a group‐level FDR threshold and a minimum delta‐PSI threshold.
#' It returns a data.frame with one row per qualifying segment.
#'
#' @param all_celltype_df A data.frame (rownames = segment IDs) (output of `test_all_groups_as()`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#'
#' @return A data.frame with one row per segment satisfying `group_fdr < fdr_thr` and `dpsi > dpsi_thr`.
#'   Columns in the returned data.frame:
#'   \describe{
#'     \item{`pv`}{Raw p‐value from the likelihood‐ratio test for the group effect (copied from `all_celltype_df$group`).}
#'     \item{`fdr`}{Benjamini‐Hochberg adjusted FDR (copied from `all_celltype_df$group_fdr`).}
#'     \item{`dpsi`}{Delta‐PSI (copied from `all_celltype_df$dpsi`).}
#'     \item{`seg_id`}{Segment ID (rownames of `all_celltype_df`).}
#'     \item{`group`}{Group label with highest mean PSI (copied from `all_celltype_df$high_state`).}
#'   }
#'   Row names of the returned data.frame are the segment IDs.
#'
#' @examples
#' \dontrun{
#' # Assume 'all_test' is output of test_all_groups_as(), with rownames = segment IDs.
#' sig_df <- select_markers_from_all_celltype_test(all_test, fdr_thr = 0.05, dpsi_thr = 0.2)
#' }
#'
#' @export
select_markers_from_all_celltype_test <- function(
    all_celltype_df,
    fdr_thr = 0.05,
    dpsi_thr = 0.1) {
  # Identify segments passing both thresholds
  keep <- (all_celltype_df$group_fdr < fdr_thr) & (all_celltype_df$dpsi > dpsi_thr)
  subset_df <- all_celltype_df[keep, , drop = FALSE]

  # Build result data.frame
  result <- data.frame(
    pv = subset_df$group,
    fdr = subset_df$group_fdr,
    dpsi = subset_df$dpsi,
    seg_id = rownames(subset_df),
    group = subset_df$high_state,
    row.names = rownames(subset_df),
    stringsAsFactors = FALSE
  )
  return(result)
}


#' Select combined marker segments from per‐group and all‐groups tests
#'
#' This function merges two sets of marker results:
#'   1. Per‐group markers (output of `select_markers`): segments deemed significant in one‐vs‐rest tests.
#'   2. All‐groups markers (output of `select_markers_from_all_celltype_test`): segments significant across all groups.
#'
#' For each group, it prioritizes segments identified by one‐vs‐rest tests; any segment not already selected but significant in the all‐groups test is labeled as a “background” marker.
#' The result contains at most `n` segments per group, with ties broken by largest |dpsi|.
#'
#' @param markers A list containing numeric matrices `pv`, `fdr`, and `dpsi`,
#'   each a matrix with rows = segments and columns = group labels (as returned by `find_marker_as()`).
#'   \describe{
#'     \item{`markers$pv`}{Matrix of raw p‐values, dimensions: ∣segments∣ × ∣groups∣.}
#'     \item{`markers$fdr`}{Matrix of Benjamini‐Hochberg adjusted FDR p‐values, same dimensions.}
#'     \item{`markers$dpsi`}{Matrix of delta‐PSI (mean PSI(group) - mean PSI(others)), same dimensions.}
#'   }
#' @param all_celltype_df A data.frame (rownames = segment IDs) (output of `test_all_groups_as()`), containing columns:
#'   \describe{
#'     \item{overdispersion}{Estimated dispersion from each segment’s GLM.}
#'     \item{group}{Raw p‐value from the likelihood‐ratio test for the `group` term.}
#'     \item{group_fdr}{Benjamini‐Hochberg adjusted FDR (across all segments).}
#'     \item{low_state}{Group label with lowest mean PSI (from `get_dpsi()`).}
#'     \item{high_state}{Group label with highest mean PSI.}
#'     \item{dpsi}{Difference in mean PSI between `high_state` and `low_state`.}
#'   }
#' @param n Integer: maximum number of markers to return per cell type (after combining) (default: `Inf` (no limit)).
#' @param fdr_thr Numeric: false discovery rate cutoff (default: 0.05).
#' @param dpsi_thr Numeric: minimum absolute delta PSI cutoff (default: 0.1).
#'
#' @return A data.frame with one row per selected segment, containing:
#'   \describe{
#'     \item{`pv`}{Raw p‐value (from either per‐group or all‐groups test).}
#'     \item{`fdr`}{Adjusted FDR (from the same test).}
#'     \item{`dpsi`}{Delta‐PSI (for the group that nominated this segment).}
#'     \item{`seg_id`}{Segment IDs (rownames).}
#'     \item{`group`}{Group label for which the segment is selected.}
#'     \item{`is_marker`}{Logical: `TRUE` if selected by the per‐group test, `FALSE` if only from the all‐groups test.}
#'   }
#'   Row names are the segment IDs. Rows are ordered within each group by descending |dpsi| and limited to `n` per group.
#'
#' @details
#' 1. Per‐group markers (priority):
#'    Call `select_markers(markers, n = Inf, fdr_thr = fdr_thr, dpsi_thr = dpsi_thr, clean_duplicates = TRUE)`
#'     to gather all segments that pass thresholds in any group. Label these as `is_marker = TRUE`.
#' 2. All‐groups markers (background):
#'    Call `select_markers_from_all_celltype_test(all_celltype_df, fdr_thr = fdr_thr, dpsi_thr = dpsi_thr)`
#'     to gather segments passing the all‐groups thresholds. Exclude any segment already in step 1. Label remaining as `is_marker = FALSE`.
#' 3. Combine and trim:
#'    For each group, among combined rows (priority first), keep at most `n` segments ranked by descending |dpsi|.
#'
#' @seealso `\link{select_markers}`, `\link{select_markers_from_all_celltype_test}`
#' @export
select_all_markers <- function(
    markers,
    all_celltype_df,
    n = Inf,
    fdr_thr = 0.05,
    dpsi_thr = 0.1) {
  # 1. Fetch all per‐group markers without limiting n
  per_group_df <- select_markers(
    markers,
    n = Inf,
    fdr_thr = fdr_thr,
    dpsi_thr = dpsi_thr,
    clean_duplicates = TRUE
  )
  # Add a column to flag these as true markers
  per_group_df$is_marker <- TRUE

  # 2. Fetch all‐groups background markers
  background_df <- select_markers_from_all_celltype_test(
    all_celltype_df,
    fdr_thr = fdr_thr,
    dpsi_thr = dpsi_thr
  )
  # Rename 'seg_id' column consistently
  # (select_markers_from_all_celltype_test already has seg_id and group)
  background_df$is_marker <- FALSE

  # 3. Exclude any segments already in per_group_df
  bg_keep <- setdiff(background_df$seg_id, per_group_df$seg_id)
  background_df <- background_df[background_df$seg_id %in% bg_keep, , drop = FALSE]

  # 4. Combine the two sets
  combined <- rbind(per_group_df, background_df)

  # 5. For each group, limit to top 'n' by |dpsi|
  # Split by group
  split_by_group <- split(combined, combined$group)
  trimmed_list <- lapply(split_by_group, function(df_group) {
    # Order by is_marker & descending |dpsi| (Prioritise select_markers())
    df_group <- df_group[order(df_group$is_marker, abs(df_group$dpsi), decreasing = TRUE), , drop = FALSE]
    # Trim to at most n rows
    if (nrow(df_group) > n) {
      df_group <- df_group[seq_len(n), , drop = FALSE]
    }
    return(df_group)
  })
  final_df <- do.call(rbind, trimmed_list)
  rownames(final_df) <- final_df$seg_id

  # 6. Re‐order rows by group (alphabetical) then |dpsi| descending
  final_df <- final_df[order(final_df$group, -abs(final_df$dpsi)), , drop = FALSE]

  return(final_df)
}


###### prepare_reference.R ######

#' Annotate coding status of segments using a GTF file
#'
#' Given a path to an Ensembl GTF (gene annotation) and an AS result list containing segment coordinates,
#'  this function determines whether each segment overlaps coding sequence (CDS) regions.
#' It assigns coding status:
#'   - 'c' if the segment is fully contained within a CDS (i.e., likely constitutive exon),
#'   - 'p' if the segment partially overlaps a CDS,
#'   - 'n' if the segment does not overlap any CDS.
#' It also flags whether the segment’s gene has at least one coding overlap.
#'
#' @param gtf_path Character: path to an Ensembl GTF file
#'   (tab‐delimited, with columns: `chr_id`, `feature` (e.g., "CDS"), `start`, `stop`, `strand`, and gene metadata).
#'   Must be readable by `read.table(..., comment.char = "#")`.
#' @param as_list SAJR::loadSAData() output list containing segment coordinates.
#' List containing at least these components:
#'   \describe{
#'     \item{`seg`}{A data.frame with columns `chr_id`, `start`, `stop`, `strand`, and `gene_id` for each segment. Row names are segment IDs.}
#'   }
#'
#' @return The input `as_list` with additional columns:
#'   \describe{
#'     \item{`seg$cod`}{Character vector of length ∣segments∣.
#'           Each element is `'c'` if the segment is completely within a CDS, `'p'` if partially overlaps, or `'n'` if no overlap.}
#'     \item{`seg$cod.gene`}{Logical vector of length ∣segments∣. `TRUE` if the segment’s gene contains at least one segment with `cod != 'n'`.}
#'   }
#'   The modified `as_list` is returned invisibly.
#'
#' @details
#' 1. Read the GTF file via `read.table(gtf_path, sep = "\t", header = FALSE, comment.char = "#")` and
#'    filter to rows where the second column (`feature`) == "CDS".
#'    Keep columns: `chr_id`, `start`, `stop`, `strand`.
#' 2. Force CDS entries to a `GRanges` using `GenomicRanges::GRanges()`.
#'    - If some `chr_id` values lack the "chr" prefix but segments use it (or vice versa), prefix or strip "chr" to match.
#' 3. Reduce the CDS ranges per chromosome/strand to merge overlapping exons (via `reduce()`).
#' 4. Build a `GRanges` for all segments from `as_list$seg`, using `chr_id`, `start`, `stop`, and `strand`.
#' 5. Use `findOverlaps(seg_gr, cds_gr, type = "any")` to find any overlap, and `findOverlaps(seg_gr, cds_gr, type = "within")` to find full containment.
#' 6. Initialise `cod = 'n'` for all segments. For indices in the "any" overlap set, set `cod = 'p'`. For indices in the "within" set, set `cod = 'c'`.
#' 7. Determine `cod.gene` by marking as `TRUE` any segment whose `gene_id` appears among those with `cod != 'n'`.
#'
#' @note
#' - Requires `GenomicRanges` for `GRanges`, `reduce()`, and `findOverlaps()`.
#' - Uses `IRanges::IRanges()` internally for range construction.
#' - If chromosome naming mismatches occur ("chr1" vs. "1"), this function attempts to harmonize by adding "chr" prefix where needed.
#'
#' @seealso `\link[GenomicRanges]{GRanges}`, `\link[GenomicRanges]{findOverlaps}`, `\link[GenomicRanges]{reduce}`
#' @export
add_is_coding_by_ens_gtf <- function(gtf_path, as_list) {
  # 1. Read GTF and keep only CDS lines
  gtf_df <- utils::read.table(gtf_path,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  # Columns: V1=chr_id, V3=feature, V4=start, V5=stop, V7=strand
  cds_df <- gtf_df[gtf_df[, 3] == "CDS", , drop = FALSE]
  colnames(cds_df)[c(1, 4, 5, 7)] <- c("chr_id", "start", "stop", "strand") # Rename columns for clarity
  cds_df <- cds_df[, c("chr_id", "start", "stop", "strand")]


  # 2. CDS Entries: GRanges
  # Harmonize chromosome naming: if some segment chr_ids have "chr" prefix but CDS lack it, or vice versa
  seg_chr <- as_list$seg$chr_id
  if (any(startsWith(seg_chr, "chr")) && !any(startsWith(cds_df$chr_id, "chr"))) {
    cds_df$chr_id <- paste0("chr", cds_df$chr_id)
  } else if (!any(startsWith(seg_chr, "chr")) && any(startsWith(cds_df$chr_id, "chr"))) {
    cds_df$chr_id <- sub("^chr", "", cds_df$chr_id)
  }

  # Build GRanges for CDS
  cds_gr <- GenomicRanges::GRanges(
    seqnames = cds_df$chr_id,
    ranges   = IRanges::IRanges(start = cds_df$start, end = cds_df$stop),
    strand   = cds_df$strand
  )


  # 3. Merge overlapping CDS ranges
  cds_gr <- GenomicRanges::reduce(cds_gr)


  # 4. Build GRanges for segments
  seg_df <- as_list$seg
  seg_gr <- GenomicRanges::GRanges(
    seqnames = seg_df$chr_id,
    ranges = IRanges::IRanges(start = seg_df$start, end = seg_df$stop),
    strand = ifelse(is.na(seg_df$strand), "*",
      ifelse(seg_df$strand == 1, "+", ifelse(seg_df$strand == -1, "-", "*"))
    )
  )


  # 5. Find overlaps: any vs. within
  ol_any <- GenomicRanges::findOverlaps(seg_gr, cds_gr, type = "any", ignore.strand = FALSE)
  ol_within <- GenomicRanges::findOverlaps(seg_gr, cds_gr, type = "within", ignore.strand = FALSE)


  # 6. Initialise coding status
  n_segs <- length(seg_gr)
  cod_stat <- rep("n", n_segs)

  # Mark partial overlaps as "p"
  ol_any_idx <- unique(S4Vectors::queryHits(ol_any))
  cod_stat[ol_any_idx] <- "p"

  # Mark full containment as "c"
  ol_within_idx <- unique(S4Vectors::queryHits(ol_within))
  cod_stat[ol_within_idx] <- "c"

  # Assign to as_list$seg$cod
  as_list$seg$cod <- cod_stat


  # 7. Determine which genes have at least one coding segment
  coding_genes <- unique(seg_df$gene_id[cod_stat != "n"])
  as_list$seg$cod.gene <- seg_df$gene_id %in% coding_genes

  return(as_list)
}

#################################


#' Aggregate single-cell counts into pseudobulk per group
#'
#' Given a `SummarizedExperiment` of single-cell assays (e.g., inclusion/exclusion counts),
#'  this function sums counts within each group (as defined by `groupby`), removes specified derivative assays,
#'  and returns a new `SummarizedExperiment` where each column corresponds to a group-level pseudobulk.
#'
#' @param se A `SummarizedExperiment` where rows are features (e.g., segments) and columns are individual cells or barcodes.
#'           Assays include `i`, `e`, `counts`, etc.
#' @param groupby Either:
#'   \itemize{
#'     \item A character vector of length `ncol(se)` giving the group label for each column, or
#'     \item One or more column names in `colData(se)`, whose values (pasted together if multiple)
#'           define the grouping factor via `get_groupby_factor()`.
#'   }
#'   After grouping, columns with the same label will be summed together.
#' @param clean_derivatives Character vector of assay names to remove before summation
#'   (e.g., `c("psi", "cpm")`). Default `c("psi", "cpm")`.
#'
#' @return A new `SummarizedExperiment` with:
#'   \describe{
#'     \item{Assays}{Each assay from `se`, except those in `clean_derivatives`, replaced by a matrix of summed values per group.}
#'     \item{rowRanges}{Copied from `rowRanges(se)`.}
#'     \item{colData}{A data.frame of group-level metadata: one row per group, including any columns from `colData(se)`
#'                      that are constant within each group, plus aggregated columns via `pseudobulk_metadata()`.}
#'   }
#'
#' @details
#' 1. Converts `groupby` into a grouping vector via `get_groupby_factor(se, groupby)`.
#' 2. Drops any assays listed in `clean_derivatives` from `se`.
#' 3. For each remaining assay, uses `visutils::calcColSums()` to sum counts for columns that share the same
#'    group label.
#' 4. Builds a metadata data.frame by calling `pseudobulk_metadata(colData(se), groupby)`, which:
#'    \enumerate{
#'      \item Splits `colData(se)` into sub-data.frames by group.
#'      \item For each group, retains columns that are constant across cells and applies aggregation functions
#'            (e.g., sum of `ncells`) for specified columns.
#'      \item Recombines group-level rows into a single data.frame, with row names matching group labels.
#'    }
#' 5. Constructs a new `SummarizedExperiment` with summed assay matrices, original `rowRanges(se)`, and
#'    the aggregated `colData` (one row per group).
#'
#' @seealso \code{\link{get_groupby_factor}}, \code{\link{pseudobulk_metadata}}, \code{\link[visutils]{calcColSums}}
#' @export
pseudobulk <- function(
    se,
    groupby,
    clean_derivatives = c("psi", "cpm")) {
  # 1. Determine grouping labels for each column
  group_factor <- get_groupby_factor(se, groupby)

  # 2. Remove specified derivative assays before summing
  to_remove <- intersect(SummarizedExperiment::assayNames(se), clean_derivatives)
  if (length(to_remove) > 0) {
    for (assay_name in to_remove) {
      SummarizedExperiment::assay(se, assay_name) <- NULL
    }
  }

  # 3. Sum each assay by group using visutils::calcColSums
  summed_assays <- list()
  for (assay_name in SummarizedExperiment::assayNames(se)) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    summed_assays[[assay_name]] <- visutils::calcColSums(mat, group_factor)
  }

  # 4. Build group-level metadata
  meta <- as.data.frame(SummarizedExperiment::colData(se))
  group_meta <- pseudobulk_metadata(meta, group_factor)
  # Reorder to match the order of columns in summed_assays (rownames of each summed matrix)
  group_meta <- group_meta[colnames(summed_assays[[1]]), , drop = FALSE]

  # 5. Construct new SummarizedExperiment
  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays = summed_assays,
    rowRanges = SummarizedExperiment::rowRanges(se),
    colData = group_meta
  )

  return(new_se)
}


#' Calculate percent spliced‐in (PSI) per segment
#'
#' Given a `SummarizedExperiment` containing inclusion (`i`) and exclusion (`e`) counts,
#'  this function computes, for each segment and sample, the percent spliced‐in: `PSI = i / (i + e)`
#' If the total counts (`i + e`) for a segment in a sample are below `min_cov`, that PSI is set to `NA`.
#'
#' @param se A `SummarizedExperiment` with assays named `i` (inclusion counts) and `e` (exclusion counts).
#'   Rows are segments; columns are pseudobulk samples.
#' @param min_cov Integer; minimum total junction coverage (`i + e`) required to compute a valid PSI.
#'   For any `i + e < min_cov`, the PSI is set to `NA`. Default is `10`.
#'
#' @return A numeric matrix of dimensions `nrow(se)` × `ncol(se)` where each entry is the PSI value
#'   for that segment and sample, or `NA` if coverage is insufficient.
#'
#' @examples
#' \dontrun{
#' # Assume 'pseudobulk_se' is a SummarizedExperiment with assays 'i' and 'e':
#' psi_mat <- calc_psi(pseudobulk_se, min_cov = 10)
#' head(psi_mat)
#' }
#'
#' @seealso \code{\link{pseudobulk}}, \code{\link{calc_cpm}}
#' @export
calc_psi <- function(se, min_cov = 10) {
  # Verify that both 'i' and 'e' assays exist
  if (!all(c("i", "e") %in% SummarizedExperiment::assayNames(se))) {
    warning("Assays 'i' and/or 'e' not found in the SummarizedExperiment; returning NULL.")
    return(NULL)
  }

  # Extract inclusion and exclusion count matrices
  inc_mat <- SummarizedExperiment::assay(se, "i")
  exc_mat <- SummarizedExperiment::assay(se, "e")

  # Compute total counts per segment × sample
  total_mat <- inc_mat + exc_mat

  # Initialize PSI matrix
  psi_mat <- inc_mat / total_mat

  # Where total < min_cov, set PSI to NA
  psi_mat[total_mat < min_cov] <- NA_real_

  return(as.matrix(psi_mat))
}


#' Calculate counts per million (CPM) per gene/feature
#'
#' Given a `SummarizedExperiment` containing raw count data (assay named `counts`),
#'   this function computes counts per million for each feature across samples: `CPM = (counts / library_size) * 1e6`
#' where `library_size` is the total sum of counts in each sample (column).
#'
#' @param se A `SummarizedExperiment` with an assay named `counts` (integer or numeric matrix).
#'   Rows are features (e.g., genes); columns are samples.
#'
#' @return A numeric matrix of the same dimensions as the `counts` assay, where each entry is the counts per million for that feature and sample.
#'   If `counts` assay is missing, returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # Assume 'se_counts' is a SummarizedExperiment with assay 'counts'
#' cpm_mat <- calc_cpm(se_counts)
#' head(cpm_mat)
#' }
#'
#' @seealso \code{\link{pseudobulk}}
#' @export
calc_cpm <- function(se) {
  # Verify that 'counts' assay exists
  if (!("counts" %in% SummarizedExperiment::assayNames(se))) {
    warning("Assay 'counts' not found in the SummarizedExperiment; returning NULL.")
    return(NULL)
  }

  # Extract raw counts matrix
  counts_mat <- SummarizedExperiment::assay(se, "counts")

  # Compute library size per sample (column sums)
  lib_sizes <- Matrix::colSums(counts_mat, na.rm = TRUE)

  # Avoid division by zero: if any lib_size == 0, set to NA to propagate NA in CPM
  lib_sizes[lib_sizes == 0] <- NA_real_

  # Compute CPM: (counts / library_size) * 1e6
  cpm_mat <- sweep(counts_mat, 2, lib_sizes, FUN = "/") * 1e6

  return(as.matrix(cpm_mat))
}


#' Aggregate column‐level metadata for pseudobulk groups
#'
#' Given a data.frame of per‐cell metadata and a grouping vector, this function aggregates metadata so that each row corresponds to one group.
#' For each group:
#'   1. Columns that are constant within the group are retained.
#'   2. For any column listed in `aggregate`, a summary function is applied (e.g., sum).
#'
#' @param meta A data.frame whose rows correspond to individual cells or samples,
#'   and whose columns contain metadata (e.g., `sample_id`, `celltype`, `ncells`, etc.).
#'   Row names should match column names of the corresponding `SummarizedExperiment` if applicable.
#' @param groupby Either:
#'   \itemize{
#'   \item A character vector of length `nrow(meta)` giving a grouping label for each row of `meta`, or
#'   \item One or more column names in `meta`. In that case, the values of those columns are pasted (row‐wise)
#'          with the default delimiter to form a grouping factor via `get_groupby_factor()`.
#'   }
#'   After grouping, columns with the same label will be summed together.
#'   See `get_groupby_factor()` for details on how the grouping vector is constructed.
#' @param aggregate A named list of functions to apply to each column within each group.
#'   For example, `list(ncells = sum, ngenes = max)` would sum the `ncells` column and take
#'    the maximum of the `ngenes` column for each group. Default is `list(ncells = sum)`.
#'
#' @return A data.frame with one row per unique group. Columns include:
#'   - Any original columns in `meta` that had the same (constant) value for all rows in that group.
#'   - For each name in `aggregate`, a new column where the corresponding function has been applied
#'     to that column’s values within each group.
#'   Row names of the returned data.frame are the group labels.
#'
#' @examples
#' \dontrun{
#' # Suppose 'cell_meta' is a data.frame of 1000 cells with columns:
#' #   sample_id, celltype, ncells, batch, library_size
#' # We want to pseudobulk by 'celltype', summing 'ncells' and taking max 'library_size':
#' group_meta <- pseudobulk_metadata(
#'   meta = cell_meta,
#'   groupby = "celltype",
#'   aggregate = list(ncells = sum, library_size = max)
#' )
#' # 'group_meta' now has one row per unique celltype.
#' }
#'
#' @seealso \code{\link{pseudobulk}}, \code{\link{get_groupby_factor}}
#' @export
pseudobulk_metadata <- function(
    meta,
    groupby,
    aggregate = list(ncells = sum)) {
  # Construct grouping factor (length = nrow(meta))
  group_factor <- get_groupby_factor(meta, groupby)

  # Split metadata by group
  meta_split <- split(meta, group_factor)

  # For each group, retain constant columns and apply aggregation functions
  aggregated_list <- lapply(meta_split, function(df_group) {
    # Determine which columns are constant within this group
    unique_counts <- sapply(df_group, function(col) length(unique(col)))
    constant_cols <- names(unique_counts)[unique_counts == 1]

    # Start result with unique values of constant columns
    result <- unique(df_group[, constant_cols, drop = FALSE])

    # Apply aggregation functions to specified columns
    for (col_name in intersect(names(aggregate), colnames(df_group))) {
      result[, col_name] <- aggregate[[col_name]](df_group[[col_name]])
    }

    return(result)
  })

  # Determine common columns across all groups (intersection of column names)
  common_cols <- Reduce(intersect, lapply(aggregated_list, colnames))

  # Combine into a single data.frame, keeping only common columns
  combined_meta <- do.call(rbind, lapply(aggregated_list, function(df) df[, common_cols, drop = FALSE]))

  return(combined_meta)
}


#' Construct a SummarizedExperiment from raw segment data and metadata
#'
#' Given a list containing segment coordinate data (`seg`) and assay matrices (`i`, `e`, etc.),
#'  this function builds a `SummarizedExperiment` by:
#' 1. Extracting `seg` (a data.frame of segment attributes) and using its columns to create
#'    `rowRanges` (`GRanges`) with genomic coordinates and metadata.
#' 2. Taking the remaining list elements as assays (each must be a matrix with rows matching `rownames(seg)`).
#' 3. Incorporating provided column metadata `col_data` as `colData`.
#'
#' @param data_list A list with at least one element named `seg`, where:
#'   - `data_list$seg` is a data.frame with columns:
#'     - `chr_id` (chromosome names matching those in the assays’ row names),
#'     - `start` (numeric start coordinate),
#'     - `stop` (numeric end coordinate),
#'     - `strand` (numeric code: `1` for '+', `-1` for '-', or `NA`/other for '*'),
#'     - plus any additional per‐segment metadata columns.
#'   - All other list elements are assay matrices (e.g., `i`, `e`, `counts`) whose row names
#'     match `rownames(data_list$seg)` and whose column names match `rownames(col_data)`.
#' @param col_data A data.frame of column‐level metadata, where each row corresponds to a column
#'   of the resulting `SummarizedExperiment`. Row names of `col_data` must match the column names
#'   of each assay matrix provided in `data_list`.
#'
#' @return A `SummarizedExperiment` with:
#'  \itemize{
#'   \item `assays`: each list element of `data_list` (except `seg`) becomes an assay, with its matrix data.
#'   \item `rowRanges`: a `GRanges` object constructed from `data_list$seg`, with ranges given by
#'     `chr_id`, `start`, `stop`, `strand` (converted to `'+'`, `'-'`, or `'*'`), and additional
#'      metadata columns stored in `elementMetadata(rowRanges)`.
#'   \item `colData`: set to the provided `col_data`.
#'  }
#'
#' @examples
#' \dontrun{
#' # Example 'raw_list' structure:
#' raw_list <- list(
#'   seg = data.frame(
#'     chr_id = c("chr1", "chr1"),
#'     start = c(100000, 200000),
#'     stop = c(100100, 200100),
#'     strand = c(1, -1),
#'     gene_id = c("GENE1", "GENE2"),
#'     row.names = c("SEG1", "SEG2"),
#'     stringsAsFactors = FALSE
#'   ),
#'   i = matrix(1:4, nrow = 2, dimnames = list(c("SEG1", "SEG2"), c("S1", "S2"))),
#'   e = matrix(5:8, nrow = 2, dimnames = list(c("SEG1", "SEG2"), c("S1", "S2")))
#' )
#' col_metadata <- data.frame(
#'   sample_id = c("S1", "S2"),
#'   celltype = c("A", "B"),
#'   row.names = c("S1", "S2"),
#'   stringsAsFactors = FALSE
#' )
#'
#' se_obj <- make_summarized_experiment(raw_list, col_metadata)
#' }
#'
#' @seealso \code{\link[GenomicRanges]{GRanges}}, \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @export
make_summarized_experiment <- function(data_list, col_data) {
  # 1. Extract segment data.frame and remove from list
  seg_df <- data_list$seg
  data_list$seg <- NULL

  # 2. Build GRanges from seg_df
  # Convert numeric strand codes to +, -, or *
  strand_vec <- ifelse(
    is.na(seg_df$strand),
    "*",
    ifelse(seg_df$strand == 1, "+",
      ifelse(seg_df$strand == -1, "-", "*")
    )
  )
  # Create GRanges; feature_id slot holds the row names of seg_df
  row_ranges <- GenomicRanges::GRanges(
    seqnames = seg_df$chr_id,
    ranges = IRanges::IRanges(start = seg_df$start, end = seg_df$stop),
    strand = strand_vec,
    feature_id = rownames(seg_df)
  )

  # 3. Remove coordinate columns from seg_df, leaving only per‐segment metadata
  meta_cols <- seg_df
  meta_cols$chr_id <- NULL
  meta_cols$start <- NULL
  meta_cols$stop <- NULL
  meta_cols$strand <- NULL

  # Attach remaining columns as elementMetadata on row_ranges
  S4Vectors::elementMetadata(row_ranges) <- meta_cols

  # 4. The remaining items in data_list are assays; ensure each is a matrix
  assay_list <- data_list
  #   (Assume user provided correct matrix dimensions and row names match seg_df)

  # 5. Construct and return SummarizedExperiment
  se_obj <- SummarizedExperiment::SummarizedExperiment(
    assays    = assay_list,
    rowRanges = row_ranges,
    colData   = col_data
  )

  return(se_obj)
}


#' Plot heatmaps of alternative splicing (PSI) and gene expression (CPM) for selected markers
#'
#' This function takes pseudobulk splicing data (`pbas`) and gene expression data (`pbge`),
#' along with a set of marker segments, and plots two side‐by‐side heatmaps:
#'   1. PSI values (percent spliced‐in) for each marker segment across groups.
#'   2. CPM values (counts per million) for corresponding genes across groups.
#'
#' @param pbas A `SummarizedExperiment` of pseudobulk splicing data. Must contain:
#'   - `assay(pbas, "i")` and `assay(pbas, "e")` (inclusion/exclusion counts).
#'   - `assay(pbas, "psi")` (percent spliced‐in) or else `psi` will be calculated internally.
#'   Rows are segment IDs; columns are group labels (matching `colData(pbas)` grouping factor).
#' @param pbge A `SummarizedExperiment` of pseudobulk gene expression data. Must contain:
#'   - `assay(pbge, "cpm")` (counts per million).
#'   Rows are gene IDs; columns are group labels matching `pbas`.
#' @param groupby A grouping label (column name in `colData(pbas)` and `colData(pbge)`)
#'   used to aggregate and order groups. Can be a single column name or a vector of column names;
#'   see `get_groupby_factor()` for details.
#' @param markers A data.frame of selected marker segments. Must contain columns:
#'   - `seg_id` (segment ID matching `rownames(pbas)`)
#'   - `dpsi` (delta‐PSI sign and magnitude; positive means inclusion in higher group)
#'   - `group` (group label for which this segment is a marker)
#'   - `is_marker` (logical; `TRUE` for top per‐group markers, `FALSE` for background markers)
#'   Optionally, `markers` may contain a column `gene_id` if plotting gene names.
#' @param psi_scale Logical; if `TRUE`, PSI values are scaled (z‐score) per segment (column of heatmap). Default is `FALSE` (raw PSI in (0,1)).
#' @param cpm_scale Logical; if `TRUE`, CPM values are standardized (z‐score) per gene (column of heatmap). Default is `TRUE`.
#' @param group_colors Optional named vector mapping group labels to colors. If `NULL`, default colors are assigned via `char2col()`.
#' @param col_palette A vector of colors (length ≥ 2) for heatmap coloring (e.g., `rev(hcl.colors(100, "RdYlBu"))`). Default is `rev(hcl.colors(100, "RdYlBu"))`.
#' @param gene_names_col Optional character. If `markers` has a column with gene names (e.g., `"name"`),
#'   specify it here so gene names appear in the column labels of the heatmap. Default is `NULL`.
#' @param ... Additional arguments passed to underlying plotting functions (`imageWithText`, `plotColorLegend2`).
#'
#' @return Invisibly returns `NULL`. The function draws two heatmaps in the current graphics device:
#'   - Left: PSI heatmap, with rows = groups, columns = marker segments.
#'   - Right: CPM heatmap, with rows = groups, columns = corresponding genes.
#'   If scaling is applied (`psi_scale` or `cpm_scale`), legends reflect z‐score ranges.
#'
#' @details
#' 1. Calls `pseudobulk()` on `pbas` and `pbge` to ensure group‐level assays (`psi` and `cpm`) exist.
#' 2. Extracts a PSI matrix (`psi_mat`) of size (n_groups × n_markers) for `markers$seg_id`.
#'    - If `markers$dpsi` is negative, flips PSI to `1 − PSI` so that all markers align directionally.
#' 3. Determines group ordering via classical Multidimensional Scaling (MDS) on `psi_mat` correlations, so that similar groups cluster together.
#' 4. Reorders `markers` by `markers$group` matching MDS order, then by `dpsi`.
#' 5. Extracts CPM matrix (`cpm_mat`) of size (n_groups × n_genes) for genes corresponding to `markers$seg_id`.
#'    - If `gene_names_col` is provided, column names become `paste0(gene_name, ":", seg_id)`.
#' 6. Scales `psi_mat` and/or `cpm_mat` if requested (`psi_scale`, `cpm_scale`).
#' 7. Prepares row annotations (`row_anns`) with:
#'    - `ct`: group label
#'    - `is_marker`: color `gray`/`white` for `TRUE`/`FALSE`
#'    - `psi_flipped`: whether a marker’s `dpsi < 0` (flipped direction; same gray/white legend).
#' 8. Draws two heatmaps side‐by‐side:
#'    - **PSI heatmap**: calls `imageWithText()` with `psi_mat`, coloring by `col_palette`, showing `ct` annotation.
#'    - **PSI legend**: if `psi_scale = FALSE`, calls `plotColorLegend2()` to display color scale (0,1) or z‐score limits.
#'    - **CPM heatmap**: calls `imageWithText()` with `cpm_mat`, similar annotations.
#'    - **CPM legend**: calls `plotColorLegend2()` to show z‐score or raw CPM scale.
#'
#' @seealso \code{\link{pseudobulk}}, \code{\link{calc_psi}}, \code{\link{calc_cpm}}, \code{\link{get_groupby_factor}}
#' @export
marker_heatmap <- function(
    pbas,
    pbge,
    groupby,
    markers,
    psi_scale = FALSE,
    cpm_scale = TRUE,
    group_colors = NULL,
    col_palette = rev(grDevices::hcl.colors(100, "RdYlBu")),
    gene_names_col = NULL,
    ...) {
  # 1. Ensure pseudobulk assays exist
  as_pb <- pseudobulk(pbas, groupby, clean_derivatives = c("psi", "cpm"))
  if (!("psi" %in% SummarizedExperiment::assayNames(as_pb))) {
    SummarizedExperiment::assay(as_pb, "psi") <- calc_psi(as_pb)
  }
  ge_pb <- pseudobulk(pbge, groupby, clean_derivatives = c())
  if ("counts" %in% SummarizedExperiment::assayNames(ge_pb)) {
    SummarizedExperiment::assay(ge_pb, "cpm") <- calc_cpm(ge_pb)
  }

  # 2. Extract PSI matrix for marker segments
  psi <- t(SummarizedExperiment::assay(as_pb, "psi")[markers$seg_id, , drop = FALSE])
  # Flip PSI for markers with negative dpsi
  psi[, markers$dpsi < 0] <- 1 - psi[, markers$dpsi < 0]

  # If is_marker is missing, default to TRUE
  if (is.null(markers$is_marker)) {
    markers$is_marker <- TRUE
  }

  # 3. Determine group ordering via MDS on PSI correlations
  cor_mat <- stats::cor(t(psi), use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  mds_coords <- stats::cmdscale(1 - cor_mat, k = 1)
  groups <- rownames(psi)[order(mds_coords[, 1])]

  # 4. Reorder markers by group (matching MDS order), then by group factor
  markers <- markers[order(markers$group), ]
  markers <- markers[order(match(markers$group, groups)), ]
  psi <- psi[groups, markers$seg_id, drop = FALSE]

  # 5. Extract CPM matrix for genes corresponding to segments
  gids <- SummarizedExperiment::rowData(as_pb)[rownames(markers), "gene_id"]
  cpm <- t(SummarizedExperiment::assay(ge_pb, "cpm")[gids, rownames(psi), drop = FALSE])
  if (!is.null(gene_names_col)) {
    gene_names <- S4Vectors::elementMetadata(SummarizedExperiment::rowRanges(ge_pb))[gids, gene_names_col]
    colnames(cpm) <- paste0(gene_names)
    colnames(psi) <- paste0(colnames(cpm), ":", colnames(psi))
  }

  # 6. Determine colorscales
  psi_zlim <- c(0, 1)
  if (psi_scale) {
    psi <- apply(psi, 2, visutils::scaleTo)
    psi_zlim <- NULL
  }
  cpm_zlim <- range(cpm, na.rm = TRUE)
  if (cpm_scale) {
    cpm <- scale(cpm)
    max_abs <- max(abs(cpm), na.rm = TRUE)
    cpm_zlim <- c(-max_abs, max_abs)
  }

  # 7. Prepare row annotations
  row_anns <- list(
    ct = markers$group,
    is_marker = as.character(markers$is_marker),
    psi_flipped = as.character(markers$dpsi < 0)
  )
  log2col <- c("TRUE" = "gray", "FALSE" = "white")
  row_ann_cols <- list(
    ct = if (is.null(group_colors)) visutils::char2col(rownames(psi)) else group_colors,
    is_marker = log2col,
    psi_flipped = log2col
  )

  # 8. Plot side-by-side with layout
  graphics::layout(matrix(c(rep(1, nrow(psi)), rep(2, nrow(psi))), ncol = 2, byrow = FALSE),
    widths = c(1, 1)
  )
  graphics::par(
    bty = "n", tcl = -0.2, mgp = c(1.3, 0.3, 0), mar = c(0, 0.5, 0, 0),
    oma = c(6, 34, 3, 1), xpd = NA
  )

  # 8a. PSI heatmap
  xaxlab <- rownames(psi)
  if (graphics::par("mfrow")[1] == 2) xaxlab <- NULL
  visutils::imageWithText(
    psi,
    "",
    colAnns = list(ct = rownames(psi)),
    rowAnns = row_anns,
    colAnnCols = list(ct = row_ann_cols$ct),
    rowAnnCols = row_ann_cols,
    xaxlab = xaxlab,
    main = "Alternative Splicing",
    col = col_palette,
    zlim = psi_zlim,
    ...
  )
  if (!psi_scale) {
    visutils::plotColorLegend2(
      graphics::grconvertX(1.02, "npc", "nfc"), 1,
      graphics::grconvertY(0.1, "npc", "nfc"),
      graphics::grconvertY(0.9, "npc", "nfc"),
      fullzlim = psi_zlim, zlim = psi_zlim,
      zfun = identity,
      z2col = function(x) visutils::num2col(x, col_palette),
      title = "PSI"
    )
  }

  # 8b. CPM heatmap
  visutils::imageWithText(
    cpm,
    "",
    colAnns = list(ct = rownames(psi)),
    rowAnns = row_anns,
    colAnnCols = list(ct = row_ann_cols$ct),
    rowAnnCols = row_ann_cols,
    main = "Gene Expression",
    col = col_palette,
    zlim = cpm_zlim,
    ...
  )
  visutils::plotColorLegend2(
    graphics::grconvertX(1.02, "npc", "nfc"), 1,
    graphics::grconvertY(0.1, "npc", "nfc"),
    graphics::grconvertY(0.9, "npc", "nfc"),
    fullzlim = cpm_zlim, zlim = cpm_zlim,
    zfun = identity,
    z2col = function(x) visutils::num2col(x, col_palette),
    title = ifelse(cpm_scale, "z-score", "CPM")
  )

  invisible(NULL)
}


###### pipeline helpers ######

#' Load intron counts as a SummarizedExperiment
#'
#' Reads a Matrix Market file (and accompanying keys) to build a `SummarizedExperiment`
#'  where rows are intron segments and columns are cell barcodes.
#' Uses `read_named_mm()` to load the counts matrix along with row and column names.
#'
#' @param path Character; file prefix for the Matrix Market data.
#'   If `path.mtx` or `path.mtx.gz` exists, `read_named_mm(path)` is used to load the sparse matrix.
#'
#' @return A `SummarizedExperiment` with:
#'   - An assay named `counts` (a sparse `dgCMatrix`) storing intron counts: rows = intron IDs, columns = barcodes.
#'   - `rowRanges`: a `GRanges` built by parsing row names (intron IDs) of the form `chr:strand:start-end`.
#'     It sets `seqnames = chr`, `ranges = IRanges(start, end)`, and `strand` as `"+"` or `"-"`.
#'   - `colData`: a data.frame with one column `barcode` (the column names of the count matrix).
#'   If no `.mtx` file is found, returns `NULL`.
#'
#' @seealso \code{\link{read_named_mm}}, \code{\link[Matrix]{readMM}}
#' @export
load_introns_as_se <- function(path) {
  # 1. Use read_named_mm() to load the sparse matrix with row/column names
  counts <- read_named_mm(path)
  if (is.null(counts)) {
    return(NULL)
  }

  # 2. Extract feature IDs and barcodes from the matrix
  feature_ids <- rownames(counts)
  barcodes <- colnames(counts)

  # 3. Parse feature IDs of the form "chr:strand:start-end"
  parts <- strsplit(feature_ids, ":")
  chr <- vapply(parts, `[`, character(1), 1)
  strand_code <- vapply(parts, `[`, character(1), 2)
  coords <- vapply(parts, `[`, character(1), 3)
  coord_split <- strsplit(coords, "-")
  start <- as.integer(vapply(coord_split, `[`, character(1), 1))
  end <- as.integer(vapply(coord_split, `[`, character(1), 2))
  strand <- ifelse(strand_code == "1", "+",
    ifelse(strand_code == "-1", "-", "*")
  )

  # 4. Build GRanges for rowRanges
  row_ranges <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end),
    strand = strand,
    feature_id = feature_ids
  )

  # 5. Construct SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays    = list(counts = counts),
    rowRanges = row_ranges,
    colData   = data.frame(barcode = barcodes, row.names = barcodes)
  )

  return(se)
}


#' Bind a list of sparse matrices by rows, aligning columns
#'
#' Given a list of matrices (possibly sparse), this function ensures that each matrix hasthe same set of column names
#'   by inserting zero‐filled sparse columns where needed, then row‐binds them into a single sparse matrix.
#'
#' @param mat_list A list of matrices or sparse matrices. Each should have identical row names
#'   but may differ in column names.
#' @param colnames_union Optional character vector of column names to enforce across all matrices.
#'   If `NULL`, defaults to the union of all column names in `mat_list`.
#'
#' @return A single sparse matrix (class `dgCMatrix`) obtained by row‐binding all matrices in
#'   `mat_list`, with columns matching `colnames_union`. Missing columns in a given matrix
#'   are filled with zero‐filled sparse columns.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' mat1 <- sparseMatrix(
#'   i = c(1, 2), j = c(1, 2), x = c(1, 2),
#'   dims = c(2, 3), dimnames = list(c("r1", "r2"), c("A", "B", "C"))
#' )
#' mat2 <- sparseMatrix(
#'   i = c(1, 2), j = c(2, 3), x = c(3, 4),
#'   dims = c(2, 3), dimnames = list(c("r1", "r2"), c("B", "C", "D"))
#' )
#' combined <- rbind_matrix(list(mat1, mat2))
#' # Result has columns A, B, C, D, with zeros filled sparsely
#' }
#'
#' @export
rbind_matrix <- function(mat_list, colnames_union = NULL) {
  # Determine union of column names if not provided
  if (is.null(colnames_union)) {
    colnames_union <- unique(unlist(lapply(mat_list, colnames)))
  }

  # Adjust each matrix to have all columns, filling missing with sparse zeros
  adjusted_list <- lapply(mat_list, function(mat) {
    # Coerce to sparse dgCMatrix if not already
    if (!inherits(mat, "dgCMatrix")) {
      mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    }
    missing_cols <- setdiff(colnames_union, colnames(mat))
    if (length(missing_cols) > 0) {
      zero_mat <- Matrix::Matrix(
        0,
        nrow = nrow(mat),
        ncol = length(missing_cols),
        sparse = TRUE,
        dimnames = list(rownames(mat), missing_cols)
      )
      mat <- cbind(mat, zero_mat)
    }
    # Reorder columns to match colnames_union
    mat[, colnames_union, drop = FALSE]
  })

  # Row-bind all adjusted sparse matrices; result remains sparse
  do.call(rbind, adjusted_list)
}


#' Subset or reorder sparse matrix columns to a specified set
#'
#' Given a matrix or sparse matrix and a target set of column names, this function returns
#' a new sparse matrix whose columns match exactly the target set. Existing columns are
#' retained (and reordered if needed), and any target column not present is added as a
#' zero‐filled sparse column.
#'
#' @param mat A matrix or sparse matrix with named columns.
#' @param target_cols Character vector specifying the desired column names (and order).
#'
#' @return A sparse matrix (`dgCMatrix`) with columns exactly `target_cols`. Any column in
#'   `target_cols` not originally present in `mat` will be inserted as a zero‐filled sparse column.
#'   Row names are unchanged.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' m <- sparseMatrix(
#'   i = c(1, 2, 1), j = c(1, 2, 3), x = c(1, 2, 3),
#'   dims = c(2, 4), dimnames = list(c("r1", "r2"), c("A", "B", "D", "E"))
#' )
#' result <- sub_cols_matrix(m, c("A", "B", "C", "D", "E"))
#' # 'C' inserted as sparse zero column; columns ordered A,B,C,D,E
#' }
#'
#' @export
sub_cols_matrix <- function(mat, target_cols) {
  # Coerce to sparse if not already
  if (!inherits(mat, "dgCMatrix")) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
  }
  missing_cols <- setdiff(target_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      sparse = TRUE,
      dimnames = list(rownames(mat), missing_cols)
    )
    mat <- cbind(mat, zero_mat)
  }
  # Reorder columns to exactly match target_cols
  mat[, target_cols, drop = FALSE]
}


#' Read a Matrix Market file with accompanying row/column keys
#'
#' Attempts to load a sparse matrix from `<prefix>.mtx` or `<prefix>.mtx.gz`,
#'  and then assigns row names from `<prefix>_key1.csv` (or `.gz`) and column names from `<prefix>_key2.csv` (or `.gz`).
#'
#' @param prefix Character; file path prefix for the Matrix Market data.
#'   If `prefix.mtx` or `prefix.mtx.gz` exists, reads it via `Matrix::readMM()`.
#'   Row names are read from `prefix_key1.csv(.gz)`, column names from `prefix_key2.csv(.gz)`.
#'
#' @return A sparse matrix (`dgCMatrix`) with row names set to feature IDs (from key1) and
#'   column names set to barcodes (from key2). If no `.mtx` file is found, returns `NULL`.
#'
#' @examples
#' \dontrun{
#' # Suppose files "data/counts.mtx" and "data/counts_key1.csv", "data/counts_key2.csv" exist:
#' mat <- read_named_mm("data/counts")
#' dim(mat)
#' head(rownames(mat))
#' head(colnames(mat))
#' }
#'
#' @seealso \code{\link[Matrix]{readMM}}, \code{\link{load_introns_as_se}}
#' @export
read_named_mm <- function(prefix) {
  # Check for .mtx file (uncompressed or gzipped)
  mtx_path <- if (file.exists(paste0(prefix, ".mtx"))) {
    paste0(prefix, ".mtx")
  } else if (file.exists(paste0(prefix, ".mtx.gz"))) {
    paste0(prefix, ".mtx.gz")
  } else {
    return(NULL)
  }

  # Read the sparse matrix
  mat <- Matrix::readMM(mtx_path)

  # Read row names from prefix_key1.csv or .gz
  key1_path <- if (file.exists(paste0(prefix, "_key1.csv"))) {
    paste0(prefix, "_key1.csv")
  } else if (file.exists(paste0(prefix, "_key1.csv.gz"))) {
    paste0(prefix, "_key1.csv.gz")
  } else {
    stop("Row key file not found: ", prefix, "_key1.csv(.gz)")
  }
  rownames(mat) <- readLines(key1_path)

  # Read column names from prefix_key2.csv or .gz
  key2_path <- if (file.exists(paste0(prefix, "_key2.csv"))) {
    paste0(prefix, "_key2.csv")
  } else if (file.exists(paste0(prefix, "_key2.csv.gz"))) {
    paste0(prefix, "_key2.csv.gz")
  } else {
    stop("Column key file not found: ", prefix, "_key2.csv(.gz)")
  }
  colnames(mat) <- readLines(key2_path)

  return(mat)
}


#' Load single-cell AS counts (inclusion/exclusion) as sparse matrices
#'
#' Given a set of segment IDs and a file prefix, this function reads inclusion (`.i.mtx`)
#'  and exclusion (`.e.mtx`) counts via `read_named_mm()`,
#'  then subsets to only the requested segments and ensures consistent columns between `i` and `e`.
#'
#' @param segs A data.frame or `SummarizedExperiment`` whose row names are the segment IDs
#'   to load. Only those segment IDs will be returned in the output matrices.
#' @param path Character; file prefix (without extension) for inclusion/exclusion data.
#'   The function expects files:
#'   - `paste0(path, ".i.mtx")` (or `.i.mtx.gz`) and accompanying keys
#'     for inclusion counts, and
#'   - `paste0(path, ".e.mtx")` (or `.e.mtx.gz`) and accompanying keys
#'     for exclusion counts.
#'   Uses `read_named_mm()`` to read each.
#'
#' @return A list with two sparse matrices of class `dgCMatrix`:
#'   - `i`: inclusion counts (rows = `rownames(segs)`, columns as in the `.i` matrix),
#'   - `e`: exclusion counts (rows = `rownames(segs)`, columns matching the inclusion matrix).
#'
#' @details
#' 1. Calls `read_named_mm(paste0(path, ".e"))`` to load the exclusion count matrix.
#' 2. Calls `read_named_mm(paste0(path, ".i"))`` to load the inclusion count matrix.
#' 3. Subsets both matrices to only rows matching `rownames(segs)``.
#' 4. Ensures that the inclusion matrix uses the same column ordering as the exclusion matrix.
#' 5. Returns a list `list(e = e_sub, i = i_sub)``.
#'
#' @seealso \code{\link{read_named_mm}}, \code{\link{load_introns_as_se}}
#' @export
load_sc_as <- function(segs, path) {
  # Read exclusion counts
  e_mat <- read_named_mm(paste0(path, ".e"))
  # Read inclusion counts
  i_mat <- read_named_mm(paste0(path, ".i"))
  # Subset to only requested segment IDs
  seg_ids <- rownames(segs)
  e_sub <- e_mat[seg_ids, , drop = FALSE]
  # Align inclusion columns to exclusion
  i_sub <- i_mat[seg_ids, colnames(e_sub), drop = FALSE]
  return(list(e = e_sub, i = i_sub))
}


#' Compute delta‐PSI (difference in percent spliced‐in) per segment
#'
#' For a `SummarizedExperiment` of splicing segments, this function first pseudobulks the data by `groupby`,
#'   computes PSI per segment (using `calc_psi()`), and then identifies for each segment
#'   the group with lowest mean PSI (`low_state`),  the group with highest mean PSI (`high_state`), and their difference (`dpsi`).
#'
#' @param data A `SummarizedExperiment` where rows are segments and columns are single‐cell or pseudobulk samples.
#'   Must contain assays `i` (inclusion) and `e` (exclusion) so that `calc_psi()` works after pseudobulking.
#' @param groupby Either:
#'   - A character vector of length `ncol(data)` giving a grouping label for each column, or
#'   - One or more column names in `colData(data)`, whose values (pasted if multiple) define the grouping factor.
#'   See `get_groupby_factor()` for details.
#' @param min_cov Integer; minimum total junction coverage (`i + e`) per pseudobulk required to compute a valid PSI.
#'                PSI values where `i + e < min_cov` are set to `NA`. Default is 50.
#'
#' @return A data.frame with one row per segment (rownames = `rownames(data)`), containing columns:
#'   - `low_state`: group label with lowest mean PSI (or `NA` if fewer than two non‐NA values),
#'   - `high_state`: group label with highest mean PSI,
#'   - `dpsi`: numeric difference `PSI[high_state] - PSI[low_state]` (or `NA` if not computable).
#'
#' @details
#' 1. The function calls `pseudobulk(data, groupby)`, which sums counts per group.
#' 2. It then computes a PSI matrix via `calc_psi()` on the pseudobulk, obtaining one PSI value per segment × group.
#' 3. For each segment (row of the PSI matrix):
#'    - Extracts non‐`NA` PSI values, sorts them by ascending PSI.
#'    - If at least two groups have non‐`NA` PSI, sets `low_state` to the name of the group with smallest PSI,
#'      `high_state` to the group with largest PSI, and `dpsi = PSI[high_state] - PSI[low_state]`.
#'    - Otherwise, leaves all three as `NA`.
#'
#' @seealso \code{\link{pseudobulk}}, \code{\link{calc_psi}}
#' @export
get_dpsi <- function(data, groupby, min_cov = 50) {
  # 1. Pseudobulk by group
  pb <- pseudobulk(data, groupby)

  # 2. Compute PSI per segment × group
  psi_mat <- calc_psi(pb, min_cov = min_cov)
  # psi_mat is a matrix (n_segments × n_groups), possibly with NA

  # 3. For each segment, find low_state, high_state, dpsi
  result_list <- apply(psi_mat, 1, function(x) {
    # x: numeric vector of PSI values for one segment across groups
    non_na <- x[!is.na(x)]
    if (length(non_na) < 2) {
      return(data.frame(
        low_state = NA_character_,
        high_state = NA_character_,
        dpsi = NA_real_
      ))
    }
    # Sort by PSI
    sorted_vals <- sort(non_na)
    low_val <- sorted_vals[1]
    high_val <- sorted_vals[length(sorted_vals)]
    # Get group names corresponding to low_val and high_val
    # 'names(non_na)' holds group labels
    low_state <- names(non_na)[which(non_na == low_val)[1]]
    high_state <- names(non_na)[which(non_na == high_val)[1]]
    data.frame(
      low_state = low_state,
      high_state = high_state,
      dpsi = high_val - low_val
    )
  })

  # 4. Combine into a data.frame with rownames = segment IDs
  result_df <- do.call(rbind, result_list)
  rownames(result_df) <- rownames(data)
  return(result_df)
}


#' Find nearest constitutive exons for a given segment
#'
#' Within a `SummarizedExperiment` of segments, this function identifies the nearest
#'  upstream and downstream “constant” exons (segments) for a given segment ID (`sid`).
#' A constant exon is defined as one where the overall PSI (summed across all samples) is `NA` or ≥ `psi_thr`.
#' The search is restricted to segments belonging to the same gene.
#'
#' @param se A `SummarizedExperiment` where `rowRanges(se)` has metadata columns:
#'   `gene_id`, `start`, `end`, `strand`, and `sites`.
#'   The assays `i` and `e` are used to compute summed PSI.
#' @param sid Character; the segment ID (rownames of `se`) for which to find nearest exons.
#' @param psi_thr Numeric between 0 and 1; constant exons are those where
#'   `sum(i) / (sum(i) + sum(e)) ≥ psi_thr` or are `NA`. Default is `0.95`.
#'
#' @return A named character vector of length 2, with names:
#'   - `up`: the segment ID of the nearest upstream constant exon (or `NA` if none),
#'   - `down`: the segment ID of the nearest downstream constant exon (or `NA` if none).
#'
#' @details
#' 1. Subset `se` to only segments belonging to the same gene as `sid`.
#' 2. Compute summed inclusion (`i_sum`) and exclusion (`e_sum`) counts across all samples:
#'     ```r
#'     i_sum <- rowSums(assay(se, "i")[same_gene, , drop = FALSE])
#'     e_sum <- rowSums(assay(se, "e")[same_gene, , drop = FALSE])
#'     psi_all <- i_sum / (i_sum + e_sum)
#'     ```
#'    Any segment with `psi_all ≥ psi_thr` or `is.na(psi_all)` is considered “constant.”
#' 3. Among constant segments, define:
#'     - “up” candidates: `sites` ∈ `c("ad", "sd")`
#'     - “down” candidates: `sites` ∈ `c("ad", "ae")`
#' 4. Order segments by genomic coordinate (`start`) in the gene’s transcriptional direction:
#'     - If `strand == "-"`, use descending `start`; otherwise ascending.
#' 5. Find the position of `sid` in that ordered list, then select the nearest constant exon
#'    before (`up`) and after (`down`) in that order.
#' 6. If the gene’s strand is negative (`"-"`), swap the `up` and `down` results in the final output.
#'
#' @seealso \code{\link{plot_segment_coverage}}, \code{\link{get_dpsi}}
#' @export
find_nearest_constant_exons <- function(se, sid, psi_thr = 0.95) {
  # 1. Verify that sid exists
  if (!(sid %in% rownames(se))) {
    stop("Segment ID '", sid, "' not found in SummarizedExperiment.")
  }

  # Extract rowRanges as a data.frame
  seg_df <- as.data.frame(SummarizedExperiment::rowRanges(se))
  # Identify the gene of interest
  gene_id <- seg_df[sid, "gene_id"]
  # Subset to all segments in this gene
  same_gene_idx <- which(seg_df$gene_id == gene_id)
  sub_df <- seg_df[same_gene_idx, , drop = FALSE]

  # 2. Compute summed inclusion/exclusion counts for these segments
  i_sub <- SummarizedExperiment::assay(se, "i")[same_gene_idx, , drop = FALSE]
  e_sub <- SummarizedExperiment::assay(se, "e")[same_gene_idx, , drop = FALSE]
  i_sum <- Matrix::rowSums(i_sub)
  e_sum <- Matrix::rowSums(e_sub)
  psi_all <- i_sum / (i_sum + e_sum)

  # 3. Identify constant segments (psi ≥ psi_thr or NA)
  constant_flag <- is.na(psi_all) | (psi_all >= psi_thr)

  # Identify up/down candidates within the gene
  sites_vec <- sub_df$sites
  up_cands_idx <- same_gene_idx[which(constant_flag & sites_vec %in% c("ad", "sd"))]
  down_cands_idx <- same_gene_idx[which(constant_flag & sites_vec %in% c("ad", "ae"))]

  # 4. Order segments by start coordinate respecting strand
  strand_val <- seg_df[sid, "strand"]
  # Get start positions for same-gene segments
  starts <- sub_df$start
  if (strand_val == "-") {
    order_rel <- order(-starts)
  } else {
    order_rel <- order(starts)
  }
  ordered_idx <- same_gene_idx[order_rel]

  # 5. Locate position of sid in the ordered list
  sid_pos <- which(ordered_idx == which(rownames(se) == sid))

  # Find nearest upstream constant exon: index in ordered_idx < sid_pos ∩ up_cands_idx
  up_id <- NA_character_
  if (length(up_cands_idx) > 0) {
    # Determine positions of up candidates in the ordered list
    up_positions <- which(ordered_idx %in% up_cands_idx)
    candidates_up <- up_positions[up_positions < sid_pos]
    if (length(candidates_up) > 0) {
      nearest_up_pos <- max(candidates_up)
      up_id <- rownames(se)[ordered_idx[nearest_up_pos]]
    }
  }

  # Find nearest downstream constant exon: index in ordered_idx > sid_pos ∩ down_cands_idx
  down_id <- NA_character_
  if (length(down_cands_idx) > 0) {
    down_positions <- which(ordered_idx %in% down_cands_idx)
    candidates_down <- down_positions[down_positions > sid_pos]
    if (length(candidates_down) > 0) {
      nearest_down_pos <- min(candidates_down)
      down_id <- rownames(se)[ordered_idx[nearest_down_pos]]
    }
  }

  # 6. If gene is on negative strand, swap up/down
  if (strand_val == "-") {
    tmp <- up_id
    up_id <- down_id
    down_id <- tmp
  }

  return(c(up = up_id, down = down_id))
}


#' Construct a contingency table from paired x and y vectors
#'
#' Given two vectors `x` and `y` of equal length, and a corresponding vector `i` of values,
#'  this function creates a matrix whose rows correspond to each unique value in `x`,
#'  columns correspond to each unique value in `y`,
#'  and entries are taken from `i` at matching `(x, y)` pairs.
#'
#' @param x A vector (atomic) of length n, representing row labels.
#' @param y A vector of length n, representing column labels.
#' @param i A vector of length n, where `i[j]` is the value to place at row = `x[j]`, column = `y[j]`.
#'
#' @return A matrix of dimension `length(unique(x)) × length(unique(y))`, with row names set to sorted unique values of `x`
#'   and column names set to sorted unique values of `y`.
#'   For each position `(r, c)`, the entry is the value from `i[j]` where `x[j] == r` and `y[j] == c`. If no such `j` exists, the entry is `NA`.
#'
#' @examples
#' \dontrun{
#' x <- c("A", "B", "A", "C")
#' y <- c("X", "Y", "X", "Z")
#' i <- c(1, 2, 3, 4)
#' # unique(x) = A, B, C; unique(y) = X, Y, Z
#' # (A, X) appears twice with values 1 and 3; the latter overwrites
#' result <- cast_xy_table(x, y, i)
#' #      X  Y  Z
#' #  A   3 NA NA
#' #  B  NA  2 NA
#' #  C  NA NA  4
#' }
#'
#' @export
cast_xy_table <- function(x, y, i) {
  # Determine sorted unique values
  ys <- sort(unique(y))
  xs <- sort(unique(x))

  # Initialize result matrix with NA
  m <- matrix(
    NA,
    nrow = length(xs),
    ncol = length(ys),
    dimnames = list(xs, ys)
  )

  # Convert to character for indexing
  x_char <- as.character(x)
  y_char <- as.character(y)

  # Fill in values
  for (j in seq_along(x_char)) {
    m[x_char[j], y_char[j]] <- i[j]
  }

  return(m)
}


#' Determine genomic plotting coordinates for a segment with flanking exons
#'
#' Given a segment ID (`sid`), a pseudobulk `SummarizedExperiment` (`pb_all`), and a gene annotation data.frame (`gene_descr`),
#'  this function finds the closest “constant” upstream and downstream exons (using `find_nearest_constant_exons()`)
#'  and returns a plotting window that extends from 50 bp upstream of the upstream exon start
#'  (or gene start if none) to 50 bp upstream of the downstream exon end (or gene end if none).
#'
#' @param sid Character scalar; the segment ID (rownames of `pb_all`) to center the plot on.
#' @param pb_all A `SummarizedExperiment` of pseudobulk splicing data. Its `rowRanges(pb_all)`
#'   must include `gene_id`, `start`, `end`, and `strand`.
#' @param gene_descr A data.frame with gene annotations, where:
#'   - Row names are gene IDs,
#'   - Columns include `start` (gene start coordinate) and `end` (gene end coordinate).
#'
#' @return A named list with two numeric elements:
#'   - `start`: the genomic coordinate to begin plotting,
#'   - `stop`: the genomic coordinate to end plotting.
#'   If no upstream exon is found, `start` is the gene start; if no downstream exon is found, `stop` is the gene end.
#'
#' @details
#' 1. Uses `find_nearest_constant_exons(pb_all, sid)` to get:
#'     - `up`: ID of the nearest upstream constant exon (or `NA`),
#'     - `down`: ID of the nearest downstream constant exon (or `NA`).
#' 2. Extracts `gene_id` for `sid` from `rowRanges(pb_all)`.
#' 3. If `up` is not `NA`, sets
#'      `start = start(rowRanges(pb_all)[up, ]) - 50`
#'    otherwise,
#'      `start = gene_descr[gene_id, "start"]`.
#' 4. If `down` is not `NA`, sets
#'      `stop = end(rowRanges(pb_all)[down, ]) - 50`
#'    otherwise,
#'      `stop = gene_descr[gene_id, "end"]`.
#' 5. Returns the two values as a named list.
#'
#' @seealso \code{\link{find_nearest_constant_exons}}, \code{\link{plot_segment_coverage}}
#' @export
get_plot_coords_for_seg <- function(sid, pb_all, gene_descr) {
  # 1. Find nearest constant exons
  nearest <- find_nearest_constant_exons(pb_all, sid)
  up_id <- nearest["up"]
  down_id <- nearest["down"]

  # 2. Extract gene_id from rowRanges
  seg <- SummarizedExperiment::rowRanges(pb_all)
  gene_id <- seg[sid, "gene_id"]

  # 3. Determine start coordinate
  if (!is.na(up_id)) {
    # use 50 bp upstream of the upstream exon's start
    up_start <- stats::start(seg[up_id])
    start_coord <- up_start - 50
  } else {
    # fallback to gene start
    start_coord <- gene_descr[gene_id, "start"]
  }

  # 4. Determine stop coordinate
  if (!is.na(down_id)) {
    # use 50 bp upstream of the downstream exon's end
    down_end <- stats::end(seg[down_id])
    stop_coord <- down_end - 50
  } else {
    # fallback to gene end
    stop_coord <- gene_descr[gene_id, "end"]
  }

  return(list(start = start_coord, stop = stop_coord))
}


#' Sum coverage and junction counts across a list of coverage objects
#'
#' Given a list of coverage objects, each containing:
#'   - `x`: a numeric vector of genomic positions,
#'   - `cov`: a numeric vector (or Rle) of coverage values at those positions,
#'   - `juncs`: a data.frame of junctions with at least three coordinate columns and a `score` column,
#' this function merges them by:
#'   1. Summing the coverage vectors (`cov`) across all objects (requires identical `x`).
#'   2. Constructing a union of all junction coordinates (first three columns of each `juncs`),
#'      initializing scores to zero, then adding each object’s `juncs$score` into the union table.
#'
#' @param cov_list A list of coverage objects. Each object must be a list with components:
#'   - `x`: a numeric vector of positions (shared across all objects),
#'   - `cov`: a numeric vector (or Rle) of coverage for those positions,
#'   - `juncs`: a data.frame with row names identifying junction IDs, at least three coordinate columns,
#'     and a numeric `score` column.
#'
#' @return A single coverage object (list) with the same structure as the inputs:
#'   - `x`: the shared positions vector,
#'   - `cov`: the summed coverage across all inputs,
#'   - `juncs`: a data.frame containing the union of all junction coordinate rows
#'     (first three columns) and a `score` column summing scores from each input.
#'
#' @details
#' 1. Single‐element shortcut
#'    If `length(cov_list) == 1`, the function simply returns `cov_list[[1]]` without modification.
#'
#' 2. Initialization
#'    Let `r` be a copy of the first coverage object, i.e. the element at index 1 of `cov_list`.
#'    Extract its `x`, `cov`, and `juncs` components for later aggregation.
#'
#' 3. Build the union of all junction coordinates
#'    - For each coverage object in `cov_list`, take the first three columns of its `juncs` data.frame.
#'    - Row‐bind them all together and then apply `unique()` to get every distinct coordinate triple.
#'    - Store these unique triples in a new data.frame called `all_juncs` and set `all_juncs$score <- 0`.
#'
#' 4. Accumulate scores from the first object
#'    - Match the row names of `r$juncs` against the row names of `all_juncs`.
#'    - For each matching row, add `r$juncs$score` into `all_juncs$score`.
#'
#' 5. Loop over remaining coverage objects
#'    For each `current` in `cov_list[2:length(cov_list)]`:
#'    - Coverage sum
#'      Check that `identical(r$x, current$x)`. If not, stop with an error. Otherwise, do
#'      ```r
#'      r$cov <- r$cov + current$cov
#'      ```
#'    - Junction score sum
#'      For each junction ID (row name) in `current$juncs`, locate the matching row in `all_juncs`
#'      (by shared row name). If found, add `current$juncs$score` to `all_juncs$score` at that row.
#'      Junctions that do not appear in `current$juncs` simply contribute zero.
#'
#' 6. Finalize
#'    After processing all objects, assign `r$juncs <- all_juncs` and return `r`.
#'
#' @export
sum_covs <- function(cov_list) {
  # If only one element, return it directly
  if (length(cov_list) == 1) {
    return(cov_list[[1]])
  }

  # 1. Initialize result as a deep copy of the first coverage object
  r <- cov_list[[1]]

  # 2. Build a data.frame of all unique junction coordinates (first three columns)
  #    across all coverage objects
  all_coords_list <- lapply(cov_list, function(cov_obj) {
    cov_obj$juncs[, 1:3, drop = FALSE]
  })
  all_juncs <- unique(do.call(rbind, all_coords_list))
  # Ensure row names of all_juncs match original junction IDs if present
  # We'll use the original row names from r$juncs for those rows that match exactly.
  # But to maintain compatibility, we'll set row names to the coordinate string.
  rownames(all_juncs) <- apply(all_juncs, 1, function(row) paste(row, collapse = ":"))

  # 3. Initialize scores to zero
  all_juncs$score <- rep(0, nrow(all_juncs))

  # 4. Add scores from the first object
  #    Map r$juncs rows into all_juncs
  #    If r$juncs row names differ from the coordinate-based rownames(all_juncs),
  #    assume rownames(r$juncs) match the coordinate paste
  matching_r <- intersect(rownames(all_juncs), rownames(r$juncs))
  if (length(matching_r) > 0) {
    all_juncs[matching_r, "score"] <-
      all_juncs[matching_r, "score"] + r$juncs[matching_r, "score"]
  }

  # 5. Loop over remaining coverage objects
  for (i in seq(2, length(cov_list))) {
    current <- cov_list[[i]]
    # Verify coverage positions match
    if (!identical(r$x, current$x)) {
      stop("Coverage positions do not match between objects.")
    }
    # Sum coverage
    r$cov <- r$cov + current$cov

    # Sum junction scores: map current$juncs rows into all_juncs
    matching_c <- intersect(rownames(all_juncs), rownames(current$juncs))
    if (length(matching_c) > 0) {
      all_juncs[matching_c, "score"] <-
        all_juncs[matching_c, "score"] + current$juncs[matching_c, "score"]
    }
  }

  # 6. Assign updated junction table back to r$juncs
  r$juncs <- all_juncs

  return(r)
}


#' Plot read coverage and splice junctions for a genomic segment across groups
#'
#' This function visualizes read coverage, junction support, PSI (percent spliced-in), and CPM for a specified segment or genomic region across multiple groups.
#' It can either accept a segment ID (`sid`) and a `SummarizedExperiment` of pseudobulk splicing (`data_as`) to compute PSI,
#'  or explicit `chr`, `start`, `stop`, and precomputed `covs` to plot raw coverage.
#' When both `data_as` and `sid` are provided, PSI boxplots are drawn;
#'  when `data_ge` is provided, CPM boxplots or plots are drawn.
#' Finally, transcript models from a `gtf` data.frame are rendered below.
#'
#' @param sid Character; segment ID (row name of `data_as`) to center on. Required if `data_as` is given.
#' @param chr Character; chromosome name (e.g., `"chr1"`). Required if not using `sid`.
#' @param start Numeric; genomic start coordinate of plot window. If `sid` is provided, can be `NULL`.
#' @param stop Numeric; genomic stop coordinate of plot window. If `sid` is provided, can be `NULL`.
#' @param covs A named list of coverage objects (each produced by `getReadCoverage()` and then `sum_covs()`),
#'   where names are group labels. If `NULL`, coverage is loaded from BAMs specified in `samples`.
#' @param celltypes Character vector of group labels (subset of those in `covs` and/or `data_as`) to plot.
#'   If `NULL` and `sid` + `data_as` provided, all groups with PSI for `sid` are used.
#' @param data_as A `SummarizedExperiment` of pseudobulk splicing data (with assays `i`, `e`, optionally `psi`).
#'   If provided along with `sid`, PSI boxplots per group will be drawn.
#' @param data_ge A `SummarizedExperiment` of pseudobulk gene‐expression data (with assay `cpm`).
#'   If provided, a CPM boxplot per group is drawn for the gene of `sid`.
#' @param groupby Character; column name in `colData(data_as)` and `colData(data_ge)` defining grouping (e.g., `"celltype"`).
#'   Required if `data_as` or `data_ge` is provided. Cannot be a length->1 vector.
#' @param barcodes A data.frame mapping `sample_id` and `barcode` to group labels (same `groupby` name). Required if raw coverage is loaded.
#' @param samples A data.frame with columns `sample_id` and `bam_path` for each sample. Required if raw coverage is loaded.
#' @param gene_descr A data.frame of gene annotations with row names = gene IDs, and columns `start`, `end`, `name`, `descr`.
#'   Required if PSI or CPM are to be labeled or if transcript trigram is needed.
#' @param scan_bam_flags A list passed to `scanBamFlags()` when loading from BAM.
#'   Default: `list(isNotPassingQualityControls = FALSE, isDuplicate = FALSE, isSupplementaryAlignment = FALSE, isSecondaryAlignment = FALSE)`.
#' @param plot_junc_only_within Logical or `NA`; if `TRUE`, only plot junctions fully within `[start, stop]`;
#'   if `FALSE`, plot junctions with at least one end in the window; if `NA`, plot all junctions. Default `NA`.
#' @param min_junc_cov_f Numeric; fraction threshold for plotting junctions by fraction of coverage. Default `0.01`.
#' @param min_junc_cov Numeric; minimum junction coverage to plot. Default `3`.
#' @param gtf A data.frame of gene/transcript annotation with columns `gene_id`, `start`, `end`, `exon.col`, `cds.col`.
#'   Only rows for the gene of `sid` are used if `sid` is provided.
#'   Must be suitable for `plotTranscripts()`.
#' @param ylim_by_junc Logical; if `TRUE`, set y‐axis limits by junction coverage only. Default `FALSE`.
#' @param ylim Numeric vector of length 2 specifying y‐axis limits for coverage plots. If `NULL`, determined automatically.
#' @param oma Numeric vector of length 4 for `par(oma = ...)`. Default `c(6, 34, 3, 1)`.
#'
#' @return Invisibly returns the updated `covs` list (with any newly loaded coverage added). Plots are drawn as a multi‐panel figure:
#'   1. Left: CPM boxplot (if `data_ge` provided).
#'   2. Middle: PSI boxplot (if `data_as` provided).
#'   3. Right: Per‐group coverage and junction plots, one row per group.
#'   4. Bottom: Transcript model plot for gene of `sid`.
#'
#' @details
#' 1. Argument validation
#'    - `groupby` must be a single column name if either `data_as` or `data_ge` is provided.
#'    - If both `start`/`stop`/`chr` are `NULL`, `sid` + `data_as` must be provided.
#'    - If raw coverage must be loaded (i.e., `covs` is `NULL`), `samples`, `barcodes`, and `groupby` are required.
#'
#' 2. PSI boxplot (when `sid` + `data_as`)
#'    - Subset `data_as` to row = `sid`.
#'    - Compute `psi <- assay(data_as, "psi")[sid, ]` if present, else call `calc_psi(data_as)[sid, ]`.
#'    - Split `psi` by group via `get_groupby_factor(data_as, groupby)` and draw horizontal boxplot via `boxplot(..., horizontal = TRUE)`.
#'
#' 3. CPM boxplot (when `data_ge` + `sid`)
#'    - Determine `gid <- rowRanges(data_as)[sid, "gene_id"]`.
#'    - Extract CPM for that gene: `cpms <- assay(data_ge, "cpm")[gid, ]`.
#'    - Split by group and draw boxplot similarly.
#'
#' 4. Coverage and junction plots
#'    - If `covs` is `NULL`, for each group:
#'      - Filter `barcodes` to those with `groupby` = group label.
#'      - Call `getReadCoverage(bam_path, chr, start, stop, strand, scanBamFlags, tagFilter = list(CB = barcodes_in_group))` for each sample,
#'        then merge via `sum_covs()`.
#'      - Store result in `covs[[group]]`.
#'    - Define layout: one column for coverage/junction per group (stacked), plus two left columns for boxplots if drawn.
#'    - For each group row:
#'      - Subset `cov <- covs[[group]]` to `[start, stop]` via `subsetCov()`.
#'      - Compute `junc_filter` based on `plot_junc_only_within`, `min_junc_cov`, and `min_junc_cov_f`.
#'      - Determine `ylim_group <- if (!is.null(ylim)) ylim else if (ylim_by_junc) c(0, max(junc_scores)) else c(0, max(cov$cov))`.
#'      - Call `plotReadCov(cov, xlim = c(start, stop), ylab = "Coverage", main = group, plot.junc.only.within = plot_junc_only_within, min.junc.cov = min_junc_cov, min.junc.cov.f = min_junc_cov_f, ylim = ylim_group, xaxt = "n")`.
#'      - Draw `abline(h = 0)`.
#'
#' 5. Transcript model
#'    - Call `plotTranscripts(gtf_subset, new = TRUE, exon.col = NA, cds.col = NA, xlim = c(start, stop))` at the bottom.
#'
#' 6. Return value
#'    Invisibly returns `covs` (a named list of per‐group coverage objects) so cached coverage can be reused.
#'
#' @seealso \code{\link{sum_covs}}, \code{\link{subset_cov}}
#' @export
plot_segment_coverage <- function(
    sid = NULL,
    chr = NULL,
    start = NULL,
    stop = NULL,
    covs = NULL,
    celltypes = NULL,
    data_as = NULL,
    data_ge = NULL,
    groupby,
    barcodes,
    samples,
    gene_descr,
    scan_bam_flags = list(
      isNotPassingQualityControls = FALSE,
      isDuplicate = FALSE,
      isSupplementaryAlignment = FALSE,
      isSecondaryAlignment = FALSE
    ),
    plot_junc_only_within = NA,
    min_junc_cov_f = 0.01,
    min_junc_cov = 3,
    gtf,
    ylim_by_junc = FALSE,
    ylim = NULL,
    oma = c(6, 34, 3, 1)) {
  # 1. Argument validation
  if ((!is.null(data_as) || !is.null(data_ge)) && length(groupby) != 1) {
    stop("`groupby` must be a single column name when using `data_as` or `data_ge`.")
  }
  if (is.null(chr) && is.null(start) && is.null(stop) && is.null(sid)) {
    stop("Either `sid` or (`chr`, `start`, `stop`) must be provided.")
  }
  if (is.null(covs) && (is.null(samples) || is.null(barcodes) || is.null(groupby))) {
    stop("To load coverage from BAM, `samples`, `barcodes`, and `groupby` are required.")
  }

  # 2. PSI preparation if sid + data_as provided
  psi <- NULL
  if (!is.null(sid) && !is.null(data_as)) {
    # Construct group factor
    seg <- as.data.frame(SummarizedExperiment::rowRanges(data_as))
    group_factor_as <- get_groupby_factor(seg, groupby)
    # Subset data_as to only this segment
    se_seg <- data_as[sid, ]
    # Compute PSI array for this segment
    if ("psi" %in% SummarizedExperiment::assayNames(se_seg)) {
      psi_vals <- SummarizedExperiment::assay(se_seg, "psi")
    } else {
      psi_vals <- calc_psi(se_seg)[1, ]
    }
    psi <- split(psi_vals, group_factor_as)
    psi <- lapply(psi, stats::na.omit)
    psi <- psi[order(sapply(psi, mean, na.rm = TRUE))]
  }

  # 3. CPM preparation if data_ge provided
  cpm <- NULL
  gid <- NULL
  if (!is.null(data_ge) && !is.null(data_as) && !is.null(sid)) {
    group_factor_ge <- get_groupby_factor(data_ge, groupby)
    gid <- seg[sid, "gene_id"]
    cpm_vals <- visutils::log10p1(SummarizedExperiment::assay(data_ge, "cpm")[gid, ])
    cpm <- split(cpm_vals, group_factor_ge)[names(psi)]
  }

  # 4. Determine genomic coordinates
  if (!is.null(sid) && !is.null(data_as)) {
    if (is.null(start)) {
      start <- gene_descr[gid, "start"]
    }
    if (is.null(stop)) {
      stop <- gene_descr[gid, "end"]
    }
    if (is.null(chr)) {
      chr <- as.character(seg[sid, "seqnames"])
    }
  }

  # 5. Prepare GTF subset for transcript model
  if (!is.null(sid)) {
    gtf <- gtf[gtf$gene_id == gid, ]
  }
  gtf$exon.col <- "black"
  gtf$cds.col <- "black"
  if (!is.null(sid)) {
    f <- gtf$start <= seg[sid, "end"] & gtf$stop >= seg[sid, "start"]
    gtf$exon.col[f] <- "red"
    gtf$cds.col[f] <- "red"
  }

  # 6. Auto‐detect celltypes if missing
  if (is.null(celltypes) && !is.null(psi)) {
    celltypes <- rev(names(psi))
  }
  if (is.null(celltypes)) {
    celltypes <- unique(barcodes[, groupby])
  }

  # 7. Build layout matrix
  #    If PSI present, allocate two left columns; if CPM, allocate one more
  n_groups <- length(celltypes)
  layout_matrix <- matrix(seq_len(1 + n_groups), ncol = 1)
  if (!is.null(psi)) {
    layout_matrix <- layout_matrix + 1
    layout_matrix <- cbind(1, layout_matrix)
  }
  if (!is.null(cpm)) {
    layout_matrix <- layout_matrix + 1
    layout_matrix <- cbind(1, layout_matrix)
  }
  if (ncol(layout_matrix) > 1) {
    layout_matrix[nrow(layout_matrix), -ncol(layout_matrix)] <- max(layout_matrix) + 1
  }
  graphics::layout(layout_matrix, widths = c(rep(1, ncol(layout_matrix) - 1), 3), heights = c(rep(1, nrow(layout_matrix) - 1), 4))
  graphics::par(bty = "n", tcl = -0.2, mgp = c(1.3, 0.3, 0), mar = c(0, 0.5, 0, 0), oma = oma, xpd = NA)

  # 8. Plot CPM boxplot if available
  if (!is.null(cpm)) {
    plot(1, t = "n", yaxs = "i", ylim = c(0.5, n_groups + 0.5), xlim = c(0, max(unlist(cpm))), yaxt = "n", xlab = "l10CPM", ylab = "")
    graphics::boxplot(cpm, horizontal = TRUE, las = 1, add = TRUE, xpd = NA, cex.axis = 2, xaxt = "n")
  }

  # 9. Plot PSI boxplot if available
  if (!is.null(psi)) {
    plot(1, t = "n", yaxs = "i", ylim = c(0.5, n_groups + 0.5), xlim = c(0, 1), yaxt = "n", xlab = "PSI", ylab = "")
    graphics::boxplot(psi, horizontal = TRUE, las = 1, add = TRUE, yaxt = if (is.null(cpm)) "s" else "n")
  }

  # 10. Coverage and junction plotting per group
  graphics::par(mar = c(0, 6, 1.1, 0), xpd = FALSE)
  for (ct in celltypes) {
    cov <- covs[[ct]]
    # Load coverage if missing or range incomplete
    if (is.null(cov) || start < cov$start || stop > cov$end) {
      cov <- list()
      for (i in seq_len(nrow(samples))) {
        sample_id <- samples$sample_id[i]
        bam_path <- samples$bam_path[i]
        tags <- barcodes$barcode[barcodes$sample_id == sample_id & !is.na(barcodes[, groupby]) & barcodes[, groupby] == ct]
        if (length(tags) == 0) next
        cov[[length(cov) + 1]] <- plotCoverage::getReadCoverage(bam_path, chr, start, stop, strand = NA, scanBamFlags = scan_bam_flags, tagFilter = list(CB = tags))
      }
      if (length(cov) > 0) {
        cov <- sum_covs(cov)
      }
      covs[[ct]] <- cov
    }
    # Subset cov to [start, stop]
    cov_sub <- subset_cov(covs[[ct]], start, stop)

    # Determine junction filter
    juncs <- cov_sub$juncs
    junc_filter <- rep(TRUE, nrow(juncs))
    if (!is.na(plot_junc_only_within) && plot_junc_only_within) {
      junc_filter <- (juncs$start >= start & juncs$end <= stop)
    } else if (!is.na(plot_junc_only_within) && !plot_junc_only_within) {
      junc_filter <- (juncs$start >= start & juncs$start <= stop) |
        (juncs$end >= start & juncs$end <= stop)
    }

    # Determine y‐axis limits
    if (is.null(ylim)) {
      raw_cov_vals <- cov_sub$cov@values
      raw_junc_vals <- juncs$score[junc_filter]
      if (ylim_by_junc) {
        y_max <- max(1, raw_junc_vals, na.rm = TRUE)
      } else {
        y_max <- max(1, raw_cov_vals, raw_junc_vals, na.rm = TRUE)
      }
      ylim_group <- c(0, y_max)
    } else {
      ylim_group <- ylim
    }

    # Plot coverage+junction
    plotCoverage::plotReadCov(cov_sub,
      xlim = c(start, stop), ylab = "Coverage", xlab = chr,
      main = ct, plot.junc.only.within = plot_junc_only_within,
      min.junc.cov = min_junc_cov, min.junc.cov.f = min_junc_cov_f,
      ylim = ylim_group, xaxt = "n"
    )
    graphics::abline(h = 0)
  }

  # 11. Transcript model plot
  graphics::par(mar = c(3, 6, 0.2, 0))
  plotCoverage::plotTranscripts(gtf, new = TRUE, exon.col = NA, cds.col = NA, xlim = c(start, stop))

  # 12. CPM vs PSI scatter if both present
  if (!is.null(psi) && !is.null(cpm)) {
    lncol <- ceiling(n_groups / 30)
    graphics::par(mar = c(3, 8 * lncol, 3, 0), xpd = NA)
    mean_cpm <- sapply(cpm, mean)
    mean_psi <- sapply(psi, mean, na.rm = TRUE)
    visutils::plotVisium(cbind(mean_cpm, mean_psi),
      ylim = c(0, 1),
      labels = names(psi), type = "p", xaxt = "s", yaxt = "s",
      pch = 16, xlab = "l10CPM", ylab = "PSI", bty = "n", cex = 2, xaxs = "r", yaxs = "r",
      legend.args = list(
        x = graphics::grconvertX(0, "ndc", "user"),
        y = graphics::grconvertY(1, "npc", "user"),
        ncol = lncol
      )
    )
  }

  # 13. Add main title across panels if sid is provided
  if (!is.null(sid)) {
    gene_info <- gene_descr[gid, ]
    graphics::mtext(paste0(sid, " ", gene_info["name"], "\n", gene_info["descr"]), side = 3, outer = TRUE)
  }

  invisible(covs)
}


#' Subset a coverage object to a specified genomic window
#'
#' Given a coverage object (a list with components `x`, `cov`, `start`, `end`, and `juncs`),
#' this function filters the coverage and junctions to only include data within `[start, stop]`.
#'
#' @param cov A coverage object (list) containing at least:
#'   - `x`: numeric vector of genomic positions,
#'   - `cov`: numeric vector or Rle of coverage values at those positions,
#'   - `start`, `end`: numeric scalars indicating the current coverage bounds,
#'   - `juncs`: data.frame of junctions with numeric columns `start` and `end`.
#' @param start Numeric; new window start coordinate.
#' @param stop Numeric; new window stop coordinate.
#'
#' @return The same coverage object `cov`, but with:
#'   - `cov$cov` and `cov$x` filtered to positions `x >= start & x <= stop`,
#'   - `cov$start` set to `start`, `cov$end` set to `stop`,
#'   - `cov$juncs` filtered to rows where `(juncs$start <= stop & juncs$end >= start)`.
#'
#' @examples
#' \dontrun{
#' # Given cov with positions 1:1000 and junctions spanning various ranges:
#' cov_window <- subset_cov(cov, start = 100, stop = 200)
#' }
#'
#' @export
subset_cov <- function(cov, start, stop) {
  # Filter coverage vector and positions
  pos_filter <- cov$x >= start & cov$x <= stop
  cov$cov <- cov$cov[pos_filter]
  cov$x <- cov$x[pos_filter]
  # Update stored start/end
  cov$start <- start
  cov$end <- stop
  # Filter junctions overlapping the window
  cov$juncs <- cov$juncs[
    cov$juncs$start <= stop & cov$juncs$end >= start, ,
    drop = FALSE
  ]
  return(cov)
}


#' Read an RDS file if it exists, otherwise return NULL
#'
#' Attempts to read an RDS file from a given path. If the file does not exist,
#' returns `NULL` without error.
#'
#' @param f Character; file path to the `.rds` file.
#'
#' @return The object returned by `readRDS(f)` if `file.exists(f) == TRUE`; otherwise `NULL`.
#'
#' @examples
#' \dontrun{
#' data_obj <- read_rds_if_exists("output/data.rds")
#' if (is.null(data_obj)) {
#'   message("No cached RDS found.")
#' }
#' }
#'
#' @export
read_rds_if_exists <- function(f) {
  if (file.exists(f)) {
    return(readRDS(f))
  }
  return(NULL)
}


#' Print an informational message with a timestamp
#'
#' Prepends `"INFO [<timestamp>]: "` to the given `text` and prints to standard output.
#'
#' @param text Character; message text or format string.
#' @param ... Additional values to append to `text` via `paste0()`.
#'
#' @return Invisibly returns `NULL`; prints `"INFO [YYYY-MM-DD HH:MM:SS]: "` followed by `text`.
#'
#' @examples
#' \dontrun{
#' log_info("Starting analysis for gene ", gene_id)
#' }
#'
#' @export
log_info <- function(text, ...) {
  msg <- paste0(text, ...)
  full_msg <- paste0("INFO [", Sys.time(), "]: ", msg)
  message(full_msg)
  invisible(NULL)
}


#' Plot MDS with group connections and labels
#'
#' Given two-dimensional coordinates for samples (e.g., output of `cmdscale()`),
#'  this function plots each point colored by its group,
#'  draws segments connecting all pairs of points within the same group,
#'  and labels each group at the mean of its points.
#'
#' @param xy A numeric matrix or data.frame with two columns (dimensions) and rows corresponding to samples.
#'   Row names should match `samples`.
#' @param celltypes A character vector giving the group label for each row of `xy`.
#'   Length must equal `nrow(xy)`.
#' @param ct2col A named vector mapping each unique group label in `celltypes` to a plotting color.
#' @param samples A vector (or factor) of sample identifiers (same order as rows of `xy`),
#'   used to choose point characters. Can be numeric or character; coerced to factor internally.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}()}.
#'
#' @return Invisibly returns `NULL`; plots a scatter of `xy` points with colored segments and labels.
#'
#' @details
#' 1. Each point `xy[i, ]` is plotted using `pch = as.numeric(factor(samples[i])) %% 25`
#'    and color `ct2col[celltypes[i]]`.
#' 2. For each distinct group `ct`, all pairs of points `xy[which(celltypes == ct), ]` are
#'    connected with line segments of color `ct2col[ct]`.
#' 3. At the mean coordinate of each group’s points, the group label `ct` is drawn.
#'
#' @examples
#' \dontrun{
#' # Suppose 'mds_coords' is a matrix of size (n_samples × 2),
#' # 'groups' is a factor of length n_samples, and 'colors' maps groups to colors.
#' plot_mds(mds_coords, groups, colors, samples = rownames(mds_coords))
#' }
#'
#' @seealso \code{\link{cmdscale}}, \code{\link{plot}}
#' @export
plot_mds <- function(xy, celltypes, ct2col, samples, ...) {
  # Validate dimensions
  if (nrow(xy) != length(celltypes) || nrow(xy) != length(samples)) {
    stop("Length of 'celltypes' and 'samples' must match nrow(xy).")
  }

  # Determine plotting character per sample (cycle through 1:24)
  pch_vals <- as.numeric(factor(samples)) %% 25
  pch_vals[pch_vals == 0] <- 25 # avoid zero, use 25 if modulo is zero

  # Initial scatter plot
  graphics::plot(xy,
    col = ct2col[celltypes],
    pch = pch_vals,
    xlab = "Dim 1",
    ylab = "Dim 2",
    ...
  )

  # Draw segments connecting all pairs within each group
  unique_ct <- unique(celltypes)
  for (ct in unique_ct) {
    idx <- which(celltypes == ct)
    if (length(idx) > 1) {
      combs <- utils::combn(idx, 2)
      graphics::segments(
        xy[combs[1, ], 1], xy[combs[1, ], 2],
        xy[combs[2, ], 1], xy[combs[2, ], 2],
        col = ct2col[ct]
      )
    }
  }

  # Label each group at its centroid
  for (ct in unique_ct) {
    idx <- which(celltypes == ct)
    centroid_x <- mean(xy[idx, 1])
    centroid_y <- mean(xy[idx, 2])
    graphics::text(centroid_x, centroid_y, labels = ct, col = ct2col[ct], cex = 1)
  }

  invisible(NULL)
}


#' Compute recommended plot dimensions for a compareCluster result
#'
#' Given a `compareClusterResult` object from `clusterProfiler`,
#'  this function returns a two‐element numeric vector giving height and width for plotting the results.
#' The height scales with the number of clusters (up to 5 rows per cluster),
#'  and the width scales with the number of clusters along the x‐axis.
#' If `ccr` is `NULL`, returns `c(2, 2)`.
#'
#' @param ccr A `compareClusterResult` object (from `clusterProfiler::compareCluster()`), or `NULL`.
#'
#' @return A numeric vector of length 2:
#'   - `height`: if `ccr` is `NULL`, `2`; otherwise `sum(pmin(5, size)) * 0.15 + 3`,
#'     where `size` is the number of terms per cluster in `ccr@compareClusterResult`.
#'   - `width`: if `ccr` is `NULL`, `2`; otherwise `length(size) * 0.6 + 7`,
#'     where `length(size)` is the number of clusters with at least one term.
#'
#' @details
#' 1. If `ccr` is `NULL`, returns `c(2, 2)`.
#' 2. Otherwise:
#'    - Extract `ccr@compareClusterResult$Cluster` to tabulate how many terms (rows) belong to each cluster.
#'    - Cap each cluster’s row count at 5 via `pmin(5, size)`.
#'    - Compute:
#'      `height = sum(pmin(5, size)) * 0.15 + 3`,
#'      `width  = length(size) * 0.6 + 7`, where `size` excludes clusters with zero terms.
#'
#' @examples
#' \dontrun{
#' dims1 <- get_dit_plot_size(ccr) # for a valid compareClusterResult
#' dims2 <- get_dit_plot_size(NULL) # returns c(2, 2)
#' }
#'
#' @seealso \code{\link{dotplot}}
#' @export
get_dit_plot_size <- function(ccr) {
  # If NULL, return default (2, 2)
  if (is.null(ccr)) {
    return(c(2, 2))
  }
  # Extract cluster assignments for each term
  clust_tbl <- table(ccr@compareClusterResult$Cluster)
  # Keep clusters with at least one term
  size <- clust_tbl[clust_tbl > 0]
  # Compute height: sum of min(5, size) * 0.15 + 3
  height <- sum(pmin(5, size)) * 0.15 + 3
  # Compute width: number of clusters * 0.6 + 7
  width <- length(size) * 0.6 + 7
  c(height = height, width = width)
}


#' Compute pairwise correlation statistics between matrix columns
#'
#' Given two numeric matrices `x` and `y` (or a single matrix `x`),
#'  this function computes Pearson correlation coefficients, sample sizes, and p-values for each pair of columns.
#'
#' @param x A numeric matrix with columns representing variables.
#' @param y Optional numeric matrix; if omitted, `y <- x`. If provided, `ncol(y)` may differ from `ncol(x)`.
#' @param ... Additional arguments passed to `stats::cor.test()` (e.g., `method`, `use`).
#'
#' @return A list with three matrices (dimensions = `ncol(x)` × `ncol(y)`):
#'   - `e`: Pearson correlation coefficients between each pair of columns,
#'   - `n`: sample sizes used in each correlation (number of non-NA pairs),
#'   - `p`: p-values from `cor.test` for each pair.
#'   Row names = `colnames(x)`, column names = `colnames(y)`.
#'
#' @details
#' 1. If `y` is missing, set `y <- x`.
#' 2. Initialize matrices `e`, `n`, and `p` with `NA` of size `ncol(x)` × `ncol(y)`.
#' 3. For each column index `i` in `x` and `j` in `y`:
#'    - Extract `x_i <- x[, i]` and `y_j <- y[, j]`.
#'    - Compute `valid <- !is.na(x_i) & !is.na(y_j)`; let `n[i, j] <- sum(valid)`.
#'    - If `n[i, j] > 2`, run
#'      ```r
#'      ct <- cor.test(x_i[valid], y_j[valid], ...)
#'      ```
#'      and set `e[i, j] <- ct$estimate`, `p[i, j] <- ct$p.value`.
#'    - Otherwise, leave `e[i, j]` and `p[i, j]` as `NA`.
#' 4. Return `list(e = e, n = n, p = p)`.
#'
#' @examples
#' \dontrun{
#' mat1 <- matrix(rnorm(100), ncol = 5)
#' mat2 <- matrix(rnorm(100), ncol = 4)
#' res <- my_cor_test(mat1, mat2)
#' dim(res$e) # 5 × 4
#' head(res$p)
#' }
#'
#' @seealso \code{\link[stats]{cor.test}}
#' @export
my_cor_test <- function(x, y = NULL, ...) {
  # If y is not provided, use x itself
  if (is.null(y)) {
    y <- x
  }

  # Ensure column names exist
  xcn <- colnames(x)
  ycn <- colnames(y)
  if (is.null(xcn)) {
    xcn <- paste0("X", seq_len(ncol(x)))
    colnames(x) <- xcn
  }
  if (is.null(ycn)) {
    ycn <- paste0("Y", seq_len(ncol(y)))
    colnames(y) <- ycn
  }

  # Initialize result matrices
  e_mat <- matrix(NA_real_,
    nrow = ncol(x), ncol = ncol(y),
    dimnames = list(xcn, ycn)
  )
  n_mat <- matrix(NA_integer_,
    nrow = ncol(x), ncol = ncol(y),
    dimnames = list(xcn, ycn)
  )
  p_mat <- matrix(NA_real_,
    nrow = ncol(x), ncol = ncol(y),
    dimnames = list(xcn, ycn)
  )

  # Loop over column pairs
  for (i in seq_len(ncol(x))) {
    for (j in seq_len(ncol(y))) {
      x_i <- x[, i]
      y_j <- y[, j]
      valid <- !is.na(x_i) & !is.na(y_j)
      n_obs <- sum(valid)
      n_mat[i, j] <- n_obs
      if (n_obs > 2) {
        ct <- stats::cor.test(x_i[valid], y_j[valid], ...)
        e_mat[i, j] <- ct$estimate
        p_mat[i, j] <- ct$p.value
      }
    }
  }

  return(list(e = e_mat, n = n_mat, p = p_mat))
}


#' Annotate exons with splice site sequences from a FASTA
#'
#' For each exon entry in a GTF data.frame, this function assigns:
#'   - `exon_number`: position of the exon in the transcript (1-based),
#'   - `transc_exon_count`: total number of exons in the transcript,
#'   - `asite`: sequence around the acceptor site,
#'   - `dsite`: sequence around the donor site.
#'
#' @param gtf A data.frame of GTF entries with columns:
#'   - `transcript_id` (identifier for each transcript),
#'   - `feature` (e.g., `"exon"`),
#'   - `chr_id` (chromosome),
#'   - `start`, `stop` (numeric genomic coordinates),
#'   - `strand` (`"+"`, `"-"`, or `"*"`).
#'   Rows may include non-exon features; only rows where `feature == "exon"` will receive exon numbers.
#' @param fa A named list of DNAString objects (FASTA) indexed by chromosome (matching `chr_id`),
#'   suitable for extracting sequence via `get_site_seq()`.
#'
#' @return The input `gtf` data.frame with added columns:
#'   - `exon_number`: integer index of the exon within its transcript (NA if not an exon),
#'   - `transc_exon_count`: total exons in that transcript (NA if not an exon),
#'   - `asite`: string sequence flanking the acceptor (5′) splice site (NA if not an exon),
#'   - `dsite`: string sequence flanking the donor (3′) splice site (NA if not an exon).
#'
#' @details
#' 1. Splits `gtf` by `transcript_id`.
#' 2. For each transcript:
#'    - Orders rows by `start` ascending (if `strand == "+"`) or descending (if `"-"`).
#'    - For rows with `feature == "exon"`, assigns `exon_number` from 1 to `n_exons` and
#'      sets `transc_exon_count = n_exons`.
#' 3. Initializes `asite` and `dsite` as `NA`.
#' 4. For exonic rows (`feature == "exon"`):
#'    - If `strand == "+"`,
#'      - `asite <- get_site_seq(chr_id, start, mars = c(12, 5), fa, rev = FALSE)`,
#'      - `dsite <- get_site_seq(chr_id, stop,  mars = c(3, 5), fa, rev = FALSE)`.
#'    - If `strand == "-"`,
#'      - `asite <- get_site_seq(chr_id, stop,  mars = c(5, 12), fa, rev = TRUE)`,
#'      - `dsite <- get_site_seq(chr_id, start, mars = c(5, 3),  fa, rev = TRUE)`.
#'
#' @seealso \code{\link{get_site_seq}}
#' @export
annotate_exons <- function(gtf, fa) {
  # Preserve original row order via index
  gtf$inx <- seq_len(nrow(gtf))

  # Split GTF by transcript_id
  split_list <- split(gtf, gtf$transcript_id)

  # Process each transcript
  split_list <- lapply(split_list, function(df) {
    # Determine exon rows
    exon_mask <- df$feature == "exon"
    if (any(exon_mask)) {
      # Order exons by start, adjusting for strand
      if (all(df$strand[exon_mask] == "-")) {
        exon_order <- order(df$start[exon_mask], decreasing = TRUE)
      } else {
        exon_order <- order(df$start[exon_mask])
      }
      n_exons <- sum(exon_mask)
      # Assign exon_number and transc_exon_count
      df$exon_number[exon_mask] <- seq_len(n_exons)[exon_order]
      df$transc_exon_count[exon_mask] <- n_exons
    }
    return(df)
  })

  # Recombine and restore original order
  gtf <- do.call(rbind, split_list)
  gtf <- gtf[order(gtf$inx), ]
  gtf$inx <- NULL

  # Initialize asite and dsite
  gtf$asite <- NA_character_
  gtf$dsite <- NA_character_

  # Compute splice site sequences for each exon row
  exon_idx <- which(gtf$feature == "exon")
  for (i in exon_idx) {
    chr <- gtf$chr_id[i]
    pos1 <- gtf$start[i]
    pos2 <- gtf$stop[i]
    strand <- gtf$strand[i]
    if (strand == "+") {
      # Acceptor: 12 bp upstream, 5 bp downstream of start
      gtf$asite[i] <- get_site_seq(chr, pos1, mars = c(12, 5), fa, rev = FALSE)
      # Donor: 3 bp upstream, 5 bp downstream of stop
      gtf$dsite[i] <- get_site_seq(chr, pos2, mars = c(3, 5), fa, rev = FALSE)
    } else if (strand == "-") {
      # Acceptor (reverse): 5 bp upstream, 12 bp downstream of stop
      gtf$asite[i] <- get_site_seq(chr, pos2, mars = c(5, 12), fa, rev = TRUE)
      # Donor (reverse): 5 bp upstream, 3 bp downstream of start
      gtf$dsite[i] <- get_site_seq(chr, pos1, mars = c(5, 3), fa, rev = TRUE)
    }
  }

  return(gtf)
}


#' Extract sequence around a splice site from a genomic FASTA
#'
#' Retrieves a nucleotide sequence flanking a specified genomic position for one or more loci.
#'
#' @param chr Character vector of chromosome names (matching names in `fa`).
#' @param pos Numeric vector of length equal to `chr`, specifying the genomic coordinate (1-based)  at the splice site for each locus.
#' @param mars Integer vector of length 2, giving the number of bases to include upstream (`mars[1]`) and downstream (`mars[2]`) of the specified position.
#' @param fa A named list of `DNAString` objects (e.g., from a BSgenome or a FASTA read into memory),
#'   where names match chromosome names in `chr`. Used to extract raw sequence.
#' @param rev Logical; if `TRUE`, the reverse complement of the extracted sequence is returned. Default is `FALSE`.
#' @param as_string Logical; if `TRUE`, returns a single concatenated string per locus; if `FALSE`,
#'   returns a list of single-character vectors. Default is `TRUE`.
#' @param to_upper Logical; if `TRUE`, converts sequences to uppercase. Default is `TRUE`.
#'
#' @return If `as_string = TRUE`, a character vector of length `length(chr)`, where each element
#'   is the concatenated sequence of length `mars[1] + mars[2] + 1` (position ± flanks).
#'   If  `as_string = FALSE`, a list of character vectors (one per locus).
#'
#' @details
#' For each index `i` in `seq_along(chr)`:
#' 1. Extract bases from `fa[[chr[i]]]` in the range `(pos[i] - mars[1]):(pos[i] + mars[2])`.
#' 2. If `rev = TRUE`, take the reverse complement of that subsequence.
#' 3. If `as_string = TRUE`, paste the bases into a single string.
#' 4. If `to_upper = TRUE`, convert to uppercase using `toupper()`.
#'
#' @examples
#' \dontrun{
#' # Suppose 'fa' is a list with FASTA sequences for chr1 and chr2:
#' seq1 <- get_site_seq("chr1", 100000, mars = c(12, 5), fa, rev = FALSE)
#' seq2 <- get_site_seq(c("chr1", "chr2"), c(100000, 200000), mars = c(5, 3), fa, rev = TRUE)
#' }
#'
#' @export
get_site_seq <- function(chr, pos, mars, fa, rev = FALSE, as_string = TRUE, to_upper = TRUE) {
  # Initialize result list
  res_list <- vector("list", length(chr))

  for (i in seq_along(chr)) {
    this_chr <- chr[i]
    this_pos <- pos[i]
    # Define start and end for extraction
    start_pos <- this_pos - mars[1]
    end_pos <- this_pos + mars[2]
    # Extract raw sequence from FASTA; assume fa[[this_chr]] is a DNAString
    subseq <- fa[[this_chr]][start_pos:end_pos]
    # If reverse complement requested, do so
    if (rev) {
      subseq <- Biostrings::reverseComplement(subseq)
    }
    # If returning as a single string
    if (as_string) {
      seq_char <- as.character(subseq)
      if (to_upper) {
        seq_char <- toupper(seq_char)
      }
      res_list[[i]] <- seq_char
    } else {
      seq_vec <- strsplit(as.character(subseq), "")[[1]]
      if (to_upper) {
        seq_vec <- toupper(seq_vec)
      }
      res_list[[i]] <- seq_vec
    }
  }

  # If as_string, collapse to character vector
  if (as_string) {
    return(unlist(res_list, use.names = FALSE))
  } else {
    return(res_list)
  }
}

#' Plot sequence logo from aligned sequences
#'
#' Generates a sequence logo from a set of aligned nucleotide sequences, using information content.
#'
#' @param seq A character vector of aligned sequences. Missing values (`NA`) are removed.
#' @param stack_height Function to compute stack height; default is `DiffLogo::informationContent`.
#' @param ... Additional arguments passed to `DiffLogo::seqLogo()`.
#'
#' @return Invisibly returns `NULL`; draws the sequence logo in the active graphics device.
#'
#' @details
#' 1. Removes any `NA` values from `seq`.
#' 2. Converts remaining sequences to uppercase.
#' 3. Computes a position weight matrix (PWM) via `DiffLogo::getPwmFromAlignment()`.
#' 4. Plots the sequence logo using `DiffLogo::seqLogo()` with the specified `stack_height`.
#'
#' @examples
#' \dontrun{
#' aligned_seqs <- c("ATGCA", "ATGGA", "ATGTA", NA, "ATGCA")
#' plot_logo(aligned_seqs)
#' }
#'
#' @seealso \code{\link[DiffLogo]{seqLogo}}, \code{\link[DiffLogo]{getPwmFromAlignment}}
#' @export
plot_logo <- function(seq, stack_height = DiffLogo::informationContent, ...) {
  # Remove NA sequences, convert to uppercase
  clean_seq <- toupper(seq[!is.na(seq)])
  # Compute PWM from alignment and plot logo
  pwm <- DiffLogo::getPwmFromAlignment(clean_seq)
  DiffLogo::seqLogo(pwm, stack_height = stack_height, ...)
}


#' Compute reverse complement of DNA sequences
#'
#' Given a character vector of nucleotide sequences, this function returns the reverse
#' complement of each sequence, preserving ambiguous bases.
#'
#' @param x A character vector of DNA sequences (e.g., "ATGC").
#'
#' @return A character vector of the same length as `x`,
#'   where each element is the reverse complement of the corresponding input sequence.
#'
#' @details
#' 1. Splits each sequence into individual characters.
#' 2. Uses `Biostrings::reverseComplement()` on a `BStringSet` constructed from `x`.
#' 3. Converts the result back to a plain character string.
#'
#' @examples
#' \dontrun{
#' seqs <- c("ATGC", "NNAGT", "ccgg")
#' rcomp(seqs)
#' # Returns c("GCAT", "ACTNN", "CCGG")
#' }
#'
#' @seealso \code{\link[Biostrings]{reverseComplement}}, \code{\link{shuffle}}
#' @export
rcomp <- function(x) {
  # Convert to uppercase for consistency
  upper_x <- toupper(x)
  # Build a BStringSet and compute reverse complement
  rc_set <- Biostrings::reverseComplement(Biostrings::BStringSet(upper_x))
  # Convert back to character vector
  as.character(rc_set)
}


#' Shuffle DNA sequences by k‐mer permutation
#'
#' Given a character vector of DNA sequences, this function generates shuffled versions
#'  of each sequence that preserve k‐mer composition, using a user‐specified k value.
#'
#' @param seqs A character vector of DNA sequences (each element may include ambiguous bases).
#' @param k Integer; length of k‐mers to preserve in the shuffle. Default is 1 (mononucleotide shuffle).
#' @param ... Additional arguments passed to `universalmotif::shuffle_sequences()`, such as `seed`.
#'
#' @return A character vector of the same length as `seqs`, where each element is a shuffled
#'    sequence preserving k‐mer composition.
#'   If `seqs[i]` contains ambiguous bases, they are shuffled along with the sequence.
#'
#' @details
#' 1. Converts each input sequence to uppercase.
#' 2. Uses `Biostrings::BStringSet()` to construct a set of DNA strings.
#' 3. Calls `universalmotif::shuffle_sequences()` on the `BStringSet`, specifying `k`.
#' 4. Converts the shuffled `BStringSet` back to a character vector.
#'
#' @examples
#' \dontrun{
#' original_seqs <- c("ATGCCGTA", "NNAGTCA", "ccggttaa")
#' shuffled <- shuffle(original_seqs, k = 2, seed = 42)
#' print(shuffled)
#' }
#'
#' @seealso \code{\link[Biostrings]{BStringSet}}, \code{\link[universalmotif]{shuffle_sequences}}
#' @export
shuffle <- function(seqs, k = 1, ...) {
  # Convert to uppercase for consistency
  upper_seqs <- toupper(seqs)
  # Build a BStringSet for shuffling
  bstrings <- Biostrings::BStringSet(upper_seqs)
  # Perform k‐mer‐preserving shuffle
  shuffled_set <- universalmotif::shuffle_sequences(bstrings, k = k, ...)
  # Convert back to character vector
  as.character(shuffled_set)
}


#' Find reverse‐complement palindrome regions in a sequence
#'
#' Identifies regions in a DNA sequence where a k‐mer and its reverse complement appear,
#' potentially extending into longer palindromic regions. Optionally shuffles the sequence
#' before searching.
#'
#' @param seq A character string representing a DNA sequence (may include ambiguous bases).
#' @param n Integer; initial k‐mer length to search for reverse complements. Default is 1.
#' @param shuffle_times Integer; number of times to apply `shuffle()` to `seq` before searching.
#'   If `shuffle_times = 0` (default), no shuffling is performed.
#' @param start Integer; 1‐based offset to add to reported coordinates (useful if `seq` is a subsequence).
#'   Default is 0.
#'
#' @return A data.frame of palindrome regions with columns:
#'   - `fpos`: forward start position of the region (1‐based, plus `start` offset),
#'   - `rpos`: reverse start position of the region (1‐based, plus `start` offset),
#'   - `len`: length of the palindrome region.
#'   Row names are constructed as `"fpos-rpos:len"`.
#'   If no palindromic regions are found, returns `NULL`.
#'
#' @details
#' 1. If `shuffle_times > 0`, replaces `seq` with its shuffled version via `shuffle(seq, k = shuffle_times)`.
#' 2. Builds a data.frame `nmers` with two columns:
#'    - `forward`: all `(length(seq) - n + 1)` k‐mers of length `n` starting at each position,
#'    - `reverse`: the reverse complement of each `forward` k‐mer.
#'    - `pos`: start positions `1:(length(seq) - n + 1)`.
#' 3. Identifies all matching pairs `(fpos, rpos)` where `forward[fpos] == reverse[rpos]`.
#' 4. Iteratively extends matches to longer palindromes (extending both ends) wherever possible.
#' 5. Constructs `rc_nmers_clean`, a data.frame of unique palindrome regions with:
#'    - `fpos`, `rpos` adjusted by `start`,
#'    - `len` = final palindrome length.
#'
#' @examples
#' \dontrun{
#' seq <- "AGCTTTCGA"
#' find_rcomp_regions(seq, n = 3)
#' # Might find the 3‐mer "AGC" at pos 1 and its reverse complement "GCT" at pos 7, etc.
#' }
#'
#' @seealso \code{\link{shuffle}}
#' @export
find_rcomp_regions <- function(seq, n = 1, shuffle_times = 0, start = 0) {
  # Optionally shuffle sequence
  if (shuffle_times > 0) {
    seq <- shuffle(seq, k = shuffle_times)
  }
  # Build all n‐mers and their reverse complements
  seq_length <- nchar(seq)
  max_pos <- seq_length - n + 1
  if (max_pos < 1) {
    return(NULL)
  }
  nmers <- data.frame(
    forward = substring(seq, 1:max_pos, n:max_pos + (n - 1)),
    pos = 1:max_pos,
    stringsAsFactors = FALSE
  )
  nmers$reverse <- rcomp(nmers$forward)

  # Find all initial matching pairs
  rc_matches <- NULL
  for (i in seq_len(nrow(nmers))) {
    matches <- which(nmers$reverse == nmers$forward[i])
    if (length(matches) > 0) {
      rc_matches <- rbind(rc_matches, data.frame(fpos = i, rpos = matches))
    }
  }
  if (is.null(rc_matches)) {
    return(NULL)
  }

  # Iteratively extend matches to longer palindromes
  rc_matches_clean <- NULL
  while (nrow(rc_matches) > 0) {
    seed <- rc_matches[1, , drop = FALSE]
    rc_matches <- rc_matches[-1, , drop = FALSE]
    fpos <- seed$fpos
    rpos <- seed$rpos
    length_k <- n
    repeat {
      next_f <- fpos + length_k
      next_r <- rpos - length_k
      if (next_f > max_pos || next_r < 1) {
        break
      }
      if (nmers$forward[next_f] == nmers$reverse[next_r]) {
        # Remove extended pair from rc_matches if present
        remove_idx <- which(rc_matches$fpos == next_f & rc_matches$rpos == next_r)
        if (length(remove_idx) > 0) {
          rc_matches <- rc_matches[-remove_idx, , drop = FALSE]
        }
        length_k <- length_k + 1
      } else {
        break
      }
    }
    rc_matches_clean <- rbind(
      rc_matches_clean,
      data.frame(
        fpos = fpos,
        rpos = rpos - length_k + n,
        len  = n + (length_k - n) * 2
      )
    )
  }

  # Filter to fpos <= rpos to avoid duplicates
  rc_matches_clean <- rc_matches_clean[rc_matches_clean$fpos <= rc_matches_clean$rpos, ]
  if (nrow(rc_matches_clean) == 0) {
    return(NULL)
  }
  # Adjust positions by offset and set row names
  rc_matches_clean$fpos <- rc_matches_clean$fpos + start
  rc_matches_clean$rpos <- rc_matches_clean$rpos + start
  rownames(rc_matches_clean) <- paste0(
    rc_matches_clean$fpos, "-", rc_matches_clean$rpos, ":", rc_matches_clean$len
  )

  return(rc_matches_clean)
}


#' Compute k-mer frequencies from a set of sequences
#'
#' Counts occurrences of all k-mers of length `k` across a vector of DNA sequences.
#'
#' @param seqs A character vector of nucleotide sequences (may include `NA`), all in uppercase.
#' @param k Integer; length of k-mers to tally. Default is 1.
#'
#' @return A named numeric vector of k-mer counts, where names are k-mer strings and values
#'   are their frequencies across all input sequences.
#'
#' @details
#' 1. Removes `NA` entries from `seqs`.
#' 2. For each sequence, extracts all substrings of length `k` starting at positions
#'    `1:(nchar(seq) - k + 1)`.
#' 3. Aggregates counts of each k-mer across all sequences and returns a named vector.
#'
#' @examples
#' \dontrun{
#' seqs <- c("ATGCA", "TGCAT", NA, "ATGCA")
#' kmer_counts <- kmer_f(seqs, k = 3)
#' # Might return counts for "ATG", "TGC", "GCA", "CAT", etc.
#' }
#'
#' @export
kmer_f <- function(seqs, k = 1) {
  seqs <- seqs[!is.na(seqs)]
  if (length(seqs) == 0) {
    return(stats::setNames(numeric(0), character(0)))
  }
  kmers <- character(0)
  for (s in seqs) {
    seq_length <- nchar(s)
    if (seq_length >= k) {
      starts <- seq(1, seq_length - k + 1, by = 1)
      kmers <- c(kmers, substring(s, starts, starts + k - 1))
    }
  }
  tbl <- table(kmers)
  stats::setNames(as.numeric(tbl), names(tbl))
}


#' Find matches or reverse-complement matches of a pattern in a sequence
#'
#' Searches for a given pattern in a DNA sequence, optionally also searching for its
#' reverse complement. Returns a data.frame of match positions and strand information.
#'
#' @param seq A single character string representing a DNA sequence.
#' @param pattern A character string or regular expression to match.
#' @param ignore_case Logical; if `TRUE`, matching is case-insensitive. Default is `TRUE`.
#' @param do_reverse Logical; if `TRUE`, also search for reverse complement of `pattern`. Default is `TRUE`.
#' @param ... Additional arguments passed to `gregexpr()`, such as `fixed = TRUE`.
#'
#' @return A data.frame with columns:
#'   - `from`: start position of each match (1-based),
#'   - `to`: end position of each match,
#'   - `dir`: `"f"` for forward matches, `"r"` for reverse complement matches.
#'   If no matches are found, returns `NULL`.
#'
#' @details
#' 1. Uses `gregexpr(pattern, seq, ignore_case = ignore_case, ...)` to find all forward matches.
#'    If `gregexpr` returns `-1`, no forward matches exist.
#' 2. If `do_reverse = TRUE`, computes the reverse complement of `pattern` via `rcomp(pattern)` and searches again for reverse matches.
#' 3. Combines forward and reverse matches into a single data.frame, sets column names, and returns.
#'
#' @examples
#' \dontrun{
#' seq <- "ATGCGTATGC"
#' hits <- find_matches(seq, "ATG")
#' print(hits)
#' # Might return positions where "ATG" (forward) or "CAT" (reverse complement) occur
#' }
#'
#' @seealso \code{\link{rcomp}}
#' @export
find_matches <- function(seq, pattern, ignore_case = TRUE, do_reverse = TRUE, ...) {
  # Forward matches
  fwd <- gregexpr(pattern, seq, ignore.case = ignore_case, ...)[[1]]
  if (fwd[1] == -1) {
    forward_df <- NULL
  } else {
    forward_df <- data.frame(
      from = fwd,
      to = fwd + attr(fwd, "match.length") - 1,
      dir = "f",
      stringsAsFactors = FALSE
    )
  }

  # Reverse complement matches
  if (do_reverse) {
    rc_pattern <- rcomp(pattern)
    rev_hits <- gregexpr(rc_pattern, seq, ignore.case = ignore_case, ...)[[1]]
    if (rev_hits[1] == -1) {
      reverse_df <- NULL
    } else {
      reverse_df <- data.frame(
        from = rev_hits,
        to = rev_hits + attr(rev_hits, "match.length") - 1,
        dir = "r",
        stringsAsFactors = FALSE
      )
    }
    combined <- rbind(forward_df, reverse_df)
    if (is.null(combined)) {
      return(NULL)
    }
    return(combined)
  }

  return(forward_df)
}
